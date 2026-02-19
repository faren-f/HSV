#!/home/bbc8731/miniconda3/envs/hsv_modulediscovery/bin/python
# this is for runing cami

#It takes protein modules, submits them to the NeDRex API to rank drugs (drug repurposing), monitors the jobs, saves predicted drugs, and validates those predictions against known drug sets.

# 1) Base directory: A folder containing multiple module subdirectories.
# Each module must contain a protein file (TSV) like cami_v3*.tsv with a gene column.
# 2) Three drug list files (text)
# Approved drugs
# Approved + experimental drugs
# Random drugs
# (One drug ID per line; used only for validation.)
# 3) NeDRex API: Online service used to build the network, rank drugs, and run validation.

import os
import json
import time
import sys
from pathlib import Path
import pandas as pd
from typing import List, Dict, Any
import requests
import logging
from logging.handlers import RotatingFileHandler
from tqdm import tqdm

# python /home/bbc8731/HSV/5_drug_repurposing/run_drug_repurposing_new.py /home/bbc8731/HSV/4_consensus_module_detection/cami_results /home/bbc8731/HSV/5_drug_repurposing/ground_truth_drugs.txt /home/bbc8731/HSV/5_drug_repurposing/randomdrugs.txt

def setup_logger(log_dir: Path) -> logging.Logger:
    """Set up a logger with both file and console handlers"""
    log_dir.mkdir(parents=True, exist_ok=True)

    logger = logging.getLogger('drug_repurposing')
    logger.setLevel(logging.INFO)

    file_handler = RotatingFileHandler(
        log_dir / 'drug_repurposing.log',
        maxBytes=10 * 1024 * 1024,
        backupCount=5
    )
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(
        logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    )

    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)  # Changed to INFO for progress messages
    console_handler.setFormatter(
        logging.Formatter('%(message)s')
    )

    logger.addHandler(file_handler)
    logger.addHandler(console_handler)
    
    return logger

class NeDRexDrugRepurposing:
    def __init__(self, base_dir: Path, true_drugs_files: Dict[str, Path], base_url: str = "https://api.nedrex.net"):
        self.base_dir = Path(base_dir)
        self.base_url = base_url
        self.logger = setup_logger(self.base_dir / 'logs')

        self.logger.debug("Initializing drug repurposing pipeline...")

        self.api_key = self.generate_api_key()
        self.headers = {'x-api-key': self.api_key}
        self.logger.debug("API key obtained successfully")

        # Load drug sets silently
        self.true_drugs_sets = {}
        for set_name, file_path in true_drugs_files.items():
            drugs = self._read_true_drugs(file_path)
            self.true_drugs_sets[set_name] = drugs
            self.logger.debug(f"Loaded {len(drugs)} drugs from {set_name}")

        self.methods = ['cami']
        self.algorithms = ['trustrank']

    def check_api_status(self):
        """Check if the API is working properly"""
        try:
            endpoint = f"{self.base_url}/licensed/list_node_collections"
            response = requests.get(endpoint, headers=self.headers)
            response.raise_for_status()
            self.logger.debug("API check successful")
            return True
        except Exception as e:
            self.logger.error(f"API check failed: {str(e)}")
            return False

    def _read_true_drugs(self, file_path: Path) -> List[str]:
        """Helper method to read and process drug IDs from file"""
        drugs = []
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                drug_id = line.split()[0]
                drug_id = drug_id.replace('drugbank.', '')
                drugs.append(drug_id)
        return drugs

    def generate_api_key(self) -> str:
        """Generate an API key by accepting the EULA"""
        endpoint = f"{self.base_url}/licensed/admin/api_key/generate"
        try:
            response = requests.post(endpoint, json={"accept_eula": True})
            response.raise_for_status()
            return response.json()
        except requests.exceptions.RequestException as e:
            self.logger.error(f"Failed to generate API key: {str(e)}")
            raise

    def check_status(self, uid: str, algorithm: str) -> Dict[str, Any]:
        """Check status with minimal output"""
        try:
            response = requests.get(
                f"{self.base_url}/licensed/{algorithm}/status",
                params={"uid": uid},
                headers=self.headers
            )
            response.raise_for_status()

            status_data = response.json()
            status = status_data.get('status', '')

            # Only log major status changes
            if status in ['completed', 'failed']:
                self.logger.info(f"Job {uid}: {status}")
                if status == 'failed':
                    self.logger.error(f"Error: {status_data.get('error', 'Unknown error')}")

            return status_data

        except Exception as e:
            self.logger.error(f"Status check failed: {str(e)}")
            return {"status": "error", "error": str(e)}

    def read_node_file(self, node_file: Path) -> List[str]:
        """Read protein IDs from a node TSV file"""
        df = pd.read_csv(node_file, sep='\t')
        return df['gene'].tolist()

    def submit_algorithm(self, algorithm: str, proteins: List[str]) -> str:
        """Submit a job with minimal output"""
        self.logger.debug(f"Submitting {algorithm} job")

        try:
            payload = {
                "seeds": proteins,
                "only_direct_drugs": True,
                "only_approved_drugs": True,
                "N": 100
            }

            endpoint = f"{self.base_url}/licensed/{algorithm}/submit"

            if algorithm == 'trustrank':
                payload["damping_factor"] = 0.85
                response = requests.post(endpoint, json=payload, headers=self.headers)
            elif algorithm == 'closeness':
                # Read graph file silently
                with open("licensed/PPDr-for-ranking.graphml", 'rb') as f:
                    files = {'graph_file': f}
                    payload_str = json.dumps({"seeds": proteins})
                    response = requests.post(
                        endpoint,
                        data={'seeds': payload_str},
                        files=files,
                        headers=self.headers
                    )

            response.raise_for_status()
            job_id = response.json()
            self.logger.debug(f"Job submitted successfully: {job_id}")
            return job_id

        except Exception as e:
            self.logger.error(f"Job submission failed: {str(e)}")
            if hasattr(response, 'text'):
                self.logger.debug(f"Error details: {response.text[:500]}")  # Only log first 500 chars
            raise

    def build_drug_protein_network(self):
        """Build and download the drug-protein network with minimal output"""
        self.logger.info("Building drug-protein network...")
        try:
            # Build network request
            endpoint = f"{self.base_url}/licensed/graph/builder"
            payload = {
                "nodes": ["drug", "protein"],
                "edges": ["drug_has_target"],
                "concise": True,
                "reviewed_proteins": [True, False],
                "drug_groups": ["approved"],
                "taxid": [9606]
            }

            # Submit build request
            response = requests.post(endpoint, json=payload, headers=self.headers)
            response.raise_for_status()
            graph_uid = response.json()
            self.graph_uid = graph_uid
            self.logger.debug(f"Graph build started with UID: {graph_uid}")

            # Monitor build status quietly
            pbar = tqdm(total=100, desc="Building network", unit="%")
            last_progress = 0
            while True:
                status_endpoint = f"{self.base_url}/licensed/graph/details/{graph_uid}"
                status = requests.get(status_endpoint, headers=self.headers).json()

                if status.get('status') == 'completed':
                    pbar.update(100 - last_progress)
                    pbar.close()
                    self.logger.info("Network build completed")

                    # Download graph silently
                    download_endpoint = f"{self.base_url}/licensed/graph/download/{graph_uid}.graphml"
                    download_response = requests.get(download_endpoint, headers=self.headers)
                    download_response.raise_for_status()

                    # Save graph without printing contents
                    graph_file = Path("licensed/PPDr-for-ranking.graphml")
                    graph_file.parent.mkdir(parents=True, exist_ok=True)
                    with open(graph_file, 'wb') as f:
                        f.write(download_response.content)
                    self.logger.info(f"Graph saved to {graph_file}")
                    return

                elif status.get('status') == 'failed':
                    pbar.close()
                    raise Exception(f"Graph building failed: {status.get('error')}")

                # Update progress bar based on status (assuming progress is available)
                current_progress = status.get('progress', {}).get('percent', last_progress)
                pbar.update(current_progress - last_progress)
                last_progress = current_progress

                # Silent wait
                time.sleep(30)

        except Exception as e:
            self.logger.error(f"Failed to build network: {str(e)}")
            raise

    def process_module(self, module_dir: Path):
        """Process module with minimal output"""
        self.logger.info(f"Processing module: {module_dir.name}")

        try:
            if not hasattr(self, 'graph_uid'):
                self.build_drug_protein_network()

            repurposing_dir = module_dir / 'drug_repurposing'
            repurposing_dir.mkdir(exist_ok=True)

            nodes_dir = module_dir

            # Process all methods silently
            for method in self.methods:
                self._process_method(method, nodes_dir, module_dir, repurposing_dir)

        except Exception as e:
            self.logger.error(f"Failed to process {module_dir.name}: {str(e)}")

    def _process_method(self, method: str, nodes_dir: Path, module_dir: Path, repurposing_dir: Path):
        """Helper method to process a single method"""

        node_file = [p for p in module_dir.iterdir() if p.is_file() and p.name.startswith("cami_v3")][0]


        

        if not node_file.exists():
            self.logger.warning(f"Node file not found for {method}: {node_file}")
            return

        self.logger.info(f"Processing method: {method}")
        proteins = self.read_node_file(node_file)
        proteins = ['uniprot.' + p if not p.startswith('uniprot.') else p for p in proteins]

        for algorithm in self.algorithms:
            self._run_algorithm(algorithm, method, proteins, repurposing_dir)

    def _run_algorithm(self, algorithm: str, method: str, proteins: List[str], repurposing_dir: Path):
        """Helper method to run a single algorithm"""
        self.logger.info(f"Running {algorithm} for {method}")

        method_dir = repurposing_dir / algorithm / method
        method_dir.mkdir(parents=True, exist_ok=True)

        try:
            uid = self.submit_algorithm(algorithm, proteins)
            self._monitor_job(uid, algorithm, method, method_dir)
        except Exception as e:
            self.logger.error(f"Error running {algorithm} for {method}: {str(e)}")
            self._save_error_info(method_dir, str(e), method, algorithm, len(proteins))

    def _run_validations(self, predicted_drugs: List[str], method_dir: Path, method: str):
        """Run validation for predicted drugs against each true drug set"""
        validation_dir = method_dir / 'validation'
        validation_dir.mkdir(parents=True, exist_ok=True)

        for drug_set_name, true_drugs in self.true_drugs_sets.items():
            set_validation_dir = validation_dir / drug_set_name
            set_validation_dir.mkdir(exist_ok=True)

            self.logger.debug(f"Running validation with {drug_set_name}")
            try:
                self._submit_validation(
                    predicted_drugs,
                    set_validation_dir,
                    method,
                    drug_set_name,
                    true_drugs
                )
            except Exception as e:
                self.logger.error(f"Validation failed for {drug_set_name}: {str(e)}")

    def _submit_validation(self, predicted_drugs: List[str], validation_dir: Path,
                           method_name: str, drug_set_name: str, true_drugs: List[str]):
        """Submit validation job with enhanced progress tracking"""
        if not predicted_drugs or not true_drugs:
            self.logger.warning(f"[{method_name}][{drug_set_name}] Skipping validation - no drugs to validate")
            return

        try:
            self.logger.info(f"[{method_name}][{drug_set_name}] Starting validation - "
                             f"Predicted: {len(predicted_drugs)}, True: {len(true_drugs)}")

            # Predicted drugs as [id, score] pairs
            predicted_pairs = [
                [f"drugbank.{d.replace('drugbank.', '')}", 1.0]
                for d in predicted_drugs
            ]

            # True drugs as simple strings
            true_drugs_clean = [
                f"drugbank.{d.split('#')[0].strip().replace('drugbank.', '')}"
                for d in true_drugs
            ]

            # Save validation input for debugging
            validation_input = {
                'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
                'method': method_name,
                'drug_set': drug_set_name,
                'num_predicted': len(predicted_pairs),
                'num_true': len(true_drugs_clean),
                'predicted_sample': predicted_pairs[:5],  # Save first 5 for reference
                'true_sample': true_drugs_clean[:5]
            }

            with open(validation_dir / 'validation_input.json', 'w') as f:
                json.dump(validation_input, f, indent=2)

            payload = {
                "test_drugs": predicted_pairs,
                "true_drugs": true_drugs_clean,
                "permutations": 10000,
                "only_approved_drugs": True
            }

            response = requests.post(
                f"{self.base_url}/licensed/validation/drug",
                json=payload,
                headers=self.headers
            )

            if response.status_code == 422:
                error_detail = response.json().get('detail', '')
                error_info = {
                    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
                    'status_code': 422,
                    'error': error_detail,
                    'validation_input': validation_input
                }
                with open(validation_dir / 'validation_error.json', 'w') as f:
                    json.dump(error_info, f, indent=2)
                self.logger.error(f"[{method_name}][{drug_set_name}] Validation API Error: {str(error_detail)[:200]}")
                return

            response.raise_for_status()
            uid = response.json()

            # Track validation progress
            validation_status = {
                'start_time': time.time(),
                'last_update': time.time(),
                'status': 'running',
                'uid': uid,
                'updates': []
            }
            
            pbar = tqdm(total=100, desc=f"Validation [{method_name}][{drug_set_name}]", unit="%")
            last_progress = 0

            while True:
                status = self._check_validation_status(uid)
                current_time = time.time()
                elapsed_minutes = (current_time - validation_status['start_time']) / 60
                current_progress = status.get('progress', {}).get('percent', last_progress)

                # Update progress bar
                pbar.update(current_progress - last_progress)
                last_progress = current_progress

                # Update every 60 seconds
                if current_time - validation_status['last_update'] >= 60:
                    status_update = {
                        'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
                        'elapsed_minutes': round(elapsed_minutes, 1),
                        'status': status.get('status', ''),
                        'progress': status.get('progress', {})
                    }
                    validation_status['updates'].append(status_update)

                    # Save progress to file
                    with open(validation_dir / 'validation_progress.json', 'w') as f:
                        json.dump(validation_status, f, indent=2)

                    self.logger.info(
                        f"[{method_name}][{drug_set_name}] "
                        f"Validation progress: {status.get('status', '')} ({elapsed_minutes:.1f} min)"
                    )
                    validation_status['last_update'] = current_time

                if status['status'] == 'completed':
                    pbar.update(100 - last_progress)
                    pbar.close()
                    results = self._get_validation_results(uid)
                    results_file = validation_dir / f"{method_name}_drug_validation.json"
                    with open(results_file, 'w') as f:
                        json.dump(results, f, indent=2)
                    self.logger.info(
                        f"[{method_name}][{drug_set_name}] "
                        f"Validation complete ({elapsed_minutes:.1f} min)"
                    )
                    break

                elif status['status'] == 'failed':
                    pbar.close()
                    error_msg = status.get('error', 'Unknown error')
                    self.logger.error(
                        f"[{method_name}][{drug_set_name}] Validation failed: {error_msg}"
                    )
                    # Save error information
                    error_info = {
                        'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
                        'elapsed_minutes': elapsed_minutes,
                        'error': error_msg,
                        'validation_status': validation_status
                    }
                    with open(validation_dir / 'validation_error.json', 'w') as f:
                        json.dump(error_info, f, indent=2)
                    break

                time.sleep(30)

        except Exception as e:
            self.logger.error(f"[{method_name}][{drug_set_name}] Validation error: {str(e)}")
            # Save exception information
            error_info = {
                'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
                'error_type': type(e).__name__,
                'error_message': str(e),
                'method': method_name,
                'drug_set': drug_set_name
            }
            with open(validation_dir / 'validation_exception.json', 'w') as f:
                json.dump(error_info, f, indent=2)

    def _check_validation_status(self, uid: str) -> Dict:
        """Check validation job status"""
        try:
            endpoint = f"{self.base_url}/licensed/validation/status"
            response = requests.get(endpoint, params={"uid": uid}, headers=self.headers)
            response.raise_for_status()
            return response.json()
        except Exception as e:
            self.logger.error(f"Error checking validation status: {str(e)}")
            return {"status": "error", "error": str(e)}

    def _get_validation_results(self, uid: str) -> Dict:
        """Get validation results"""
        endpoint = f"{self.base_url}/licensed/validation/status"
        response = requests.get(endpoint, params={"uid": uid}, headers=self.headers)
        response.raise_for_status()
        return response.json()

    def _monitor_job(self, uid: str, algorithm: str, method: str, method_dir: Path):
        """Monitor job progress with status updates"""
        start_time = time.time()
        last_update = 0
        pbar = tqdm(total=100, desc=f"{algorithm} for {method}", unit="%")

        while True:
            status_response = self.check_status(uid, algorithm)
            status = status_response.get('status', '')
            current_progress = status_response.get('progress', {}).get('percent', 0)

            current_time = time.time()
            elapsed_minutes = (current_time - start_time) / 60

            # Update progress bar
            pbar.update(current_progress - pbar.n)

            # Update every 60 seconds
            if current_time - last_update >= 60:
                self.logger.info(f"{algorithm} for {method}: {status} ({elapsed_minutes:.1f} min)")
                last_update = current_time

            if status == 'completed':
                pbar.update(100 - pbar.n)
                pbar.close()
                self._handle_completed_job(status_response, method_dir, method)
                break
            elif status == 'failed':
                pbar.close()
                self.logger.error(f"Job failed: {status_response.get('error', 'Unknown error')}")
                break

            time.sleep(30)

    def _handle_completed_job(self, status_response: Dict, method_dir: Path, method: str):
        """Helper method to handle completed job"""
        if 'results' in status_response:
            results = status_response['results']
            results_file = method_dir / "results.json"
            with open(results_file, 'w') as f:
                json.dump(results, f, indent=2)

            predicted_drugs = self._extract_predicted_drugs(results)
            self.logger.info(f"Got {len(predicted_drugs)} predicted drugs")

            self._run_validations(predicted_drugs, method_dir, method)

    def _extract_predicted_drugs(self, results: Dict) -> List[str]:
        """Helper method to extract predicted drugs from results"""
        if isinstance(results, dict) and 'drugs' in results:
            return [drug['drug_name'] for drug in results['drugs']]
        elif isinstance(results, list):
            return [drug['drug_name'] for drug in results
                    if isinstance(drug, dict) and 'drug_name' in drug]
        return []

    def _save_error_info(self, method_dir: Path, error: str, method: str,
                         algorithm: str, num_proteins: int):
        """Helper method to save error information"""
        error_info = {
            'error': error,
            'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
            'method': method,
            'algorithm': algorithm,
            'num_proteins': num_proteins
        }
        with open(method_dir / "error.json", 'w') as f:
            json.dump(error_info, f, indent=2)


# gets base_dir + 3 drug-list files from command line.
def main():
    if len(sys.argv) != 4:  # Updated to accept one more argument
        print("Usage: script.py <base_directory> <approved_drugs_file> <random_drugs_file>")
        sys.exit(1)

    base_dir = Path(sys.argv[1])
    true_drugs_files = {
        "approved_drugs": Path(sys.argv[2]),
        "random_drugs": Path(sys.argv[3])  # Added the random drugs file
    }

    try:
        pipeline = NeDRexDrugRepurposing(base_dir, true_drugs_files) # Create pipeline object: sets up logging, gets an API key, and loads the “true drug” lists.

        if not pipeline.check_api_status():
            pipeline.logger.error("API is not available")
            sys.exit(1)

        module_dirs = [d for d in base_dir.glob("*") if d.is_dir()]

        for module_dir in tqdm(module_dirs, desc="Overall progress"):
            pipeline.process_module(module_dir)

    except Exception as e:
        logging.error(f"Fatal error: {str(e)}", exc_info=True)
        sys.exit(1)

if __name__ == "__main__":
    main()

