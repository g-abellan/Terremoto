import subprocess
import sys
import os

def run_script(script_name):
    print(f"\n==================================================")
    print(f"Running: {script_name}")
    print(f"==================================================")
    result = subprocess.run([sys.executable, script_name], capture_output=False, text=True)
    if result.returncode != 0:
        print(f"Error: {script_name} failed with exit code {result.returncode}")
        sys.exit(result.returncode)
    print(f"Success: {script_name} completed.")

def main():
    # Ensure directories exist
    os.makedirs("../Reports", exist_ok=True)
    
    # Change working directory to the Sources directory to preserve relative paths
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    
    # Execution pipeline in order
    pipeline = [
        "fetch_data.py",
        "etas_simulation.py",
        "geometry_analysis.py",
        "field_theory_analysis.py",
        "coupling_analysis.py",
        "run_forecast.py"
    ]
    
    for script in pipeline:
        run_script(script)
        
    print("\n==================================================")
    print("Pipeline Execution Completed Successfully!")
    print("All datasets and figures have been updated.")
    print("==================================================")

if __name__ == "__main__":
    main()
