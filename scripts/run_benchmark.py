#!/usr/bin/env python3
import os
import sys
import subprocess
import glob
import re
import csv
import argparse
import matplotlib.pyplot as plt
from pathlib import Path

# --- Configuration (Based on User Inputs) ---
DEFAULT_BENCHMARK_DIR = "/public_new/project/ACH20250801_genome_T2T/test_pseudo/software/benchmark"
DEFAULT_GENOME = "/public_new/project/ACH20250801_genome_T2T/test_pseudo/whale_pseudogene_multithreading/01.genome/dolphin_chr.fa"
DEFAULT_EXECUTABLE = "/public_new/project/ACH20250801_genome_T2T/test_pseudo/software/EasyPseudogene_v1/bin/easypseudogene"
DEFAULT_WISE_CONFIG = "/public/software/wise2.4.1/wisecfg/"
DEFAULT_WISE_PATH = "/public/software/wise2.4.1/src/bin/"
DEFAULT_THREADS = 96

def run_command(cmd):
    """Run a command and capture stderr which contains the 'time' output."""
    print(f"Executing: {' '.join(cmd)}")
    
    # We use /usr/bin/time -v, which prints to stderr
    result = subprocess.run(
        cmd, 
        capture_output=True, 
        text=True
    )
    
    if result.returncode != 0:
        print(f"Error executing command: {result.stderr}")
        return None
    
    return result.stderr

def parse_time_output(stderr_output):
    """Extract Elapsed Time and Max RSS from GNU time output."""
    time_data = {}
    
    # Elapsed (wall clock) time (h:mm:ss or m:ss)
    time_match = re.search(r'Elapsed \(wall clock\) time \(h:mm:ss or m:ss\):\s+([0-9:.]+)', stderr_output)
    if time_match:
        time_str = time_match.group(1)
        parts = time_str.split(':')
        seconds = 0
        if len(parts) == 3: # h:mm:ss
            seconds = int(parts[0])*3600 + int(parts[1])*60 + float(parts[2])
        elif len(parts) == 2: # m:ss
            seconds = int(parts[0])*60 + float(parts[1])
        time_data['seconds'] = seconds
        
    # Maximum resident set size (kbytes)
    mem_match = re.search(r'Maximum resident set size \(kbytes\):\s+(\d+)', stderr_output)
    if mem_match:
        # Convert kbytes to GB
        time_data['memory_gb'] = int(mem_match.group(1)) / 1024 / 1024
        
    return time_data

def get_seq_count(filename):
    """Extract sequence count from filename (e.g., proteins_1000.fa -> 1000)."""
    match = re.search(r'proteins_(\d+)\.fa', filename)
    if match:
        return int(match.group(1))
    return 0

def plot_benchmark(csv_file):
    """Generate plots from the CSV results."""
    data = []
    with open(csv_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            data.append({
                'count': int(row['seq_count']),
                'seconds': float(row['seconds']),
                'memory': float(row['memory_gb'])
            })
    
    # Sort by sequence count
    data.sort(key=lambda x: x['count'])
    
    counts = [x['count'] for x in data]
    times = [x['seconds'] for x in data]
    mems = [x['memory'] for x in data]
    
    # Setup styles
    plt.style.use('ggplot')
    
    # 1. Runtime Plot
    plt.figure(figsize=(10, 6))
    plt.plot(counts, times, 'o-', color='#e74c3c', linewidth=2, markersize=8)
    plt.title('EasyPseudogene Scalability: Runtime vs Input Size')
    plt.xlabel('Number of Input Protein Sequences')
    plt.ylabel('Runtime (seconds)')
    plt.grid(True)
    
    # If range is large, use log scale
    if max(counts) / (min(counts) + 1) > 100:
        plt.xscale('log')
        plt.yscale('log')
        plt.title('EasyPseudogene Scalability (Log-Log Scale)')
        
    plt.savefig('benchmark_runtime.png')
    print("Generated benchmark_runtime.png")
    
    # 2. Memory Plot
    plt.figure(figsize=(10, 6))
    plt.bar([str(c) for c in counts], mems, color='#3498db', alpha=0.7)
    plt.title('EasyPseudogene Memory Usage')
    plt.xlabel('Number of Input Protein Sequences')
    plt.ylabel('Peak Memory Usage (GB)')
    plt.grid(axis='y')
    plt.savefig('benchmark_memory.png')
    print("Generated benchmark_memory.png")

def main():
    parser = argparse.ArgumentParser(description="Run EasyPseudogene benchmark")
    parser.add_argument("--input-dir", default=DEFAULT_BENCHMARK_DIR, help="Directory with input fasta files")
    parser.add_argument("--skip-run", action="store_true", help="Skip running and just plot existing CSV")
    args = parser.parse_args()
    
    results_csv = "benchmark_results.csv"
    
    if not args.skip_run:
        # Find input files
        pattern = os.path.join(args.input_dir, "proteins_*.fa")
        files = glob.glob(pattern)
        
        if not files:
            print(f"No files found matching {pattern}")
            sys.exit(1)
            
        print(f"Found {len(files)} input files to test.")
        
        results = []
        
        # Open CSV for writing
        with open(results_csv, 'w', newline='') as csvfile:
            fieldnames = ['seq_count', 'seconds', 'memory_gb', 'filename']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            
            for prot_file in sorted(files, key=get_seq_count):
                seq_count = get_seq_count(os.path.basename(prot_file))
                print(f"\n--- Processing {os.path.basename(prot_file)} ({seq_count} seqs) ---")
                
                output_dir = f"bench_out_{seq_count}"
                
                cmd = [
                    "/usr/bin/time", "-v",
                    DEFAULT_EXECUTABLE,
                    "--proteins", prot_file,
                    "--genome", DEFAULT_GENOME,
                    "--threads", str(DEFAULT_THREADS),
                    "--wiseconfig", DEFAULT_WISE_CONFIG,
                    "--wisepath", DEFAULT_WISE_PATH,
                    "--output", output_dir
                ]
                
                stderr = run_command(cmd)
                
                if stderr:
                    metrics = parse_time_output(stderr)
                    if metrics:
                        row = {
                            'seq_count': seq_count,
                            'seconds': metrics.get('seconds', 0),
                            'memory_gb': metrics.get('memory_gb', 0),
                            'filename': os.path.basename(prot_file)
                        }
                        writer.writerow(row)
                        csvfile.flush() # Save progress
                        print(f"Completed in {row['seconds']:.2f}s, Peak Mem: {row['memory_gb']:.2f}GB")
                    else:
                        print("Failed to parse time output")
    
    if os.path.exists(results_csv):
        print("\nUsing results from benchmark_results.csv to generate plots...")
        plot_benchmark(results_csv)
    else:
        print("No results CSV found.")

if __name__ == "__main__":
    main()
