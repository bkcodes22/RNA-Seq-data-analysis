from flask import Flask, render_template, request, redirect, url_for, flash, Response, jsonify, send_from_directory
import sys
import os
import subprocess
import re
import threading
import json
import time
from werkzeug.utils import secure_filename

app = Flask(__name__)
app.secret_key = 'supersecretkey'

# Global storage for processing status and terminal output
processing_status = {}
terminal_output = {}

# Fixed STAR and MultiQC executable paths
STAR_EXECUTABLE = "/home/compbio001/anaconda3/envs/STAR/bin/STAR-avx2"
MULTIQC_EXECUTABLE = "/home/compbio001/anaconda3/envs/multiqc/bin/multiqc"

UPLOAD_FOLDER = 'uploads'
os.makedirs(UPLOAD_FOLDER, exist_ok=True)

# Store for job IDs
job_counter = 0
job_lock = threading.Lock()

def find_multiqc_html(output_dir):
    """Find MultiQC HTML file in output directory"""
    for root, dirs, files in os.walk(output_dir):
        for file in files:
            if file.endswith('.html') and 'multiqc' in file.lower():
                return os.path.join(root, file)
    return None

def process_files_async(job_id, genome_dir, output_dir, threads, fastq_pairs):
    """Process files asynchronously and send updates"""
    global processing_status, terminal_output
    
    processing_status[job_id] = {
        'status': 'running',
        'current_file': 0,
        'total_files': len(fastq_pairs),
        'messages': [],
        'multiqc_html': None
    }
    terminal_output[job_id] = []
    
    try:
        total_files = len(fastq_pairs)
        
        for idx, (forward_path, reverse_path) in enumerate(fastq_pairs, 1):
            # Check if files are compressed
            is_compressed = forward_path.endswith('.gz') or (reverse_path and reverse_path.endswith('.gz'))
            
            if reverse_path:
                # Paired-end
                forward_name = os.path.basename(forward_path)
                reverse_name = os.path.basename(reverse_path)
                filename = f"{forward_name} + {reverse_name}"
                
                processing_status[job_id]['messages'].append(
                    f"🔄 [{idx}/{total_files}] Processing paired-end: {forward_name} + {reverse_name}"
                )
                processing_status[job_id]['current_file'] = idx
                
                # Extract sample name from forward read
                sample_name = os.path.basename(forward_path)
                for pattern in [r'_R1.*', r'_1.*', r'\.1\..*', r'_forward.*', r'_f\..*']:
                    sample_name = re.sub(pattern, '', sample_name, flags=re.IGNORECASE)
                if sample_name.endswith('.gz'):
                    sample_name = sample_name[:-3]
                sample_name = os.path.splitext(sample_name)[0]
                output_prefix = os.path.join(output_dir, sample_name)

                star_command = [
                    STAR_EXECUTABLE,
                    "--runThreadN", threads,
                    "--genomeDir", genome_dir,
                    "--readFilesIn", forward_path, reverse_path,
                    "--outFileNamePrefix", output_prefix,
                    "--outSAMtype", "BAM", "SortedByCoordinate",
                    "--quantMode", "GeneCounts"
                ]
                
                if is_compressed:
                    star_command.extend(["--readFilesCommand", "zcat"])
            else:
                # Single-end
                filename = os.path.basename(forward_path)
                
                processing_status[job_id]['messages'].append(
                    f"🔄 [{idx}/{total_files}] Processing single-end: {filename}"
                )
                processing_status[job_id]['current_file'] = idx
                
                sample_name = os.path.basename(forward_path)
                if sample_name.endswith('.gz'):
                    sample_name = sample_name[:-3]
                sample_name = os.path.splitext(sample_name)[0]
                output_prefix = os.path.join(output_dir, sample_name)

                star_command = [
                    STAR_EXECUTABLE,
                    "--runThreadN", threads,
                    "--genomeDir", genome_dir,
                    "--readFilesIn", forward_path,
                    "--outFileNamePrefix", output_prefix,
                    "--outSAMtype", "BAM", "SortedByCoordinate",
                    "--quantMode", "GeneCounts"
                ]
                
                if is_compressed:
                    star_command.extend(["--readFilesCommand", "zcat"])

            # Run command and capture output
            try:
                process = subprocess.Popen(
                    star_command,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    universal_newlines=True,
                    bufsize=1
                )
                
                for line in iter(process.stdout.readline, ''):
                    if line:
                        terminal_output[job_id].append(line.rstrip())
                        # Keep only last 1000 lines
                        if len(terminal_output[job_id]) > 1000:
                            terminal_output[job_id].pop(0)
                
                process.wait()
                
                if process.returncode == 0:
                    processing_status[job_id]['messages'].append(
                        f"✅ [{idx}/{total_files}] Completed: {filename}"
                    )
                else:
                    processing_status[job_id]['messages'].append(
                        f"❌ [{idx}/{total_files}] Error processing {filename}: exit code {process.returncode}"
                    )
            except Exception as e:
                processing_status[job_id]['messages'].append(
                    f"❌ [{idx}/{total_files}] Error processing {filename}: {str(e)}"
                )

        # Run MultiQC
        processing_status[job_id]['messages'].append("📊 Running MultiQC...")
        try:
            process = subprocess.Popen(
                [MULTIQC_EXECUTABLE, output_dir, "-o", output_dir],
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                universal_newlines=True,
                bufsize=1
            )
            
            for line in iter(process.stdout.readline, ''):
                if line:
                    terminal_output[job_id].append(line.rstrip())
                    if len(terminal_output[job_id]) > 1000:
                        terminal_output[job_id].pop(0)
            
            process.wait()
            
            if process.returncode == 0:
                # Find MultiQC HTML file
                multiqc_html = find_multiqc_html(output_dir)
                if multiqc_html:
                    processing_status[job_id]['multiqc_html'] = multiqc_html
                    processing_status[job_id]['messages'].append(
                        f"📊 MultiQC analysis completed. Report available."
                    )
                else:
                    processing_status[job_id]['messages'].append("📊 MultiQC analysis completed.")
            else:
                processing_status[job_id]['messages'].append(
                    f"❌ Error running MultiQC: exit code {process.returncode}"
                )
        except Exception as e:
            processing_status[job_id]['messages'].append(f"❌ Error running MultiQC: {str(e)}")
        
        processing_status[job_id]['status'] = 'completed'
        
    except Exception as e:
        processing_status[job_id]['status'] = 'error'
        processing_status[job_id]['messages'].append(f"❌ Fatal error: {str(e)}")

def pair_fastq_files(fastq_files):
    """
    Pair forward and reverse FASTQ files based on common naming patterns.
    Returns a list of tuples: [(forward_path, reverse_path)] for paired-end,
    or [(single_path, None)] for single-end files that couldn't be paired.
    """
    paired_files = []
    unpaired_files = list(fastq_files)
    
    # Common patterns for paired-end naming
    patterns = [
        (r'_R1(\.|_)', r'_R2(\.|_)'),  # sample_R1.fastq, sample_R2.fastq
        (r'_1(\.|_)', r'_2(\.|_)'),    # sample_1.fastq, sample_2.fastq
        (r'\.1\.', r'\.2\.'),          # sample.1.fastq, sample.2.fastq
        (r'_forward(\.|_)', r'_reverse(\.|_)'),  # sample_forward.fastq, sample_reverse.fastq
        (r'_f(\.|_)', r'_r(\.|_)'),    # sample_f.fastq, sample_r.fastq
        (r'\.f\.', r'\.r\.'),          # sample.f.fastq, sample.r.fastq
    ]
    
    # Try to pair files
    while unpaired_files:
        file1 = unpaired_files.pop(0)
        basename1 = os.path.basename(file1)
        matched = False
        
        for pattern_fwd, pattern_rev in patterns:
            if re.search(pattern_fwd, basename1, re.IGNORECASE):
                # Found a forward read, look for matching reverse
                potential_reverse = re.sub(pattern_fwd, pattern_rev, basename1, flags=re.IGNORECASE)
                
                for file2 in unpaired_files:
                    basename2 = os.path.basename(file2)
                    if basename2 == potential_reverse or re.search(pattern_rev, basename2, re.IGNORECASE):
                        # Check if they share the same base name
                        base1 = re.sub(pattern_fwd, '', basename1, flags=re.IGNORECASE)
                        base2 = re.sub(pattern_rev, '', basename2, flags=re.IGNORECASE)
                        if base1 == base2:
                            paired_files.append((file1, file2))
                            unpaired_files.remove(file2)
                            matched = True
                            break
                
                if matched:
                    break
        
        if not matched:
            # Couldn't find a pair, treat as single-end
            paired_files.append((file1, None))
    
    return paired_files

@app.route("/", methods=["GET", "POST"])
def index():
    if request.method == "POST":
        genome_dir = request.form.get("genome_dir").strip()
        output_dir = request.form.get("output_dir").strip()
        threads = request.form.get("threads", "8").strip()
        mode = request.form.get("mode")
        read_type = request.form.get("read_type", "single")

        if not os.path.exists(genome_dir):
            return jsonify({"error": "Invalid genome directory path!"}), 400

        os.makedirs(output_dir, exist_ok=True)

        fastq_pairs = []

        if mode == "single":
            if read_type == "single":
                uploaded_file = request.files.get("fastq_file")
                if uploaded_file and (uploaded_file.filename.endswith(".fastq") or uploaded_file.filename.endswith(".fastq.gz")):
                    filename = secure_filename(uploaded_file.filename)
                    file_path = os.path.join(UPLOAD_FOLDER, filename)
                    uploaded_file.save(file_path)
                    fastq_pairs = [(file_path, None)]
                else:
                    return jsonify({"error": "Please upload a valid FASTQ file (.fastq or .fastq.gz)."}), 400
            else:  # paired-end
                uploaded_file_1 = request.files.get("fastq_file_1")
                uploaded_file_2 = request.files.get("fastq_file_2")
                valid_extensions = (".fastq", ".fastq.gz")
                if uploaded_file_1 and any(uploaded_file_1.filename.endswith(ext) for ext in valid_extensions) and \
                   uploaded_file_2 and any(uploaded_file_2.filename.endswith(ext) for ext in valid_extensions):
                    filename_1 = secure_filename(uploaded_file_1.filename)
                    filename_2 = secure_filename(uploaded_file_2.filename)
                    file_path_1 = os.path.join(UPLOAD_FOLDER, filename_1)
                    file_path_2 = os.path.join(UPLOAD_FOLDER, filename_2)
                    uploaded_file_1.save(file_path_1)
                    uploaded_file_2.save(file_path_2)
                    fastq_pairs = [(file_path_1, file_path_2)]
                else:
                    return jsonify({"error": "Please upload both forward and reverse FASTQ files (.fastq or .fastq.gz)."}), 400

        elif mode == "batch":
            input_dir = request.form.get("input_dir").strip()
            if not os.path.exists(input_dir):
                return jsonify({"error": "Invalid input directory path!"}), 400
            
            fastq_files = [os.path.join(input_dir, f) for f in os.listdir(input_dir) 
                          if f.endswith(".fastq") or f.endswith(".fastq.gz")]
            
            if not fastq_files:
                return jsonify({"error": "No FASTQ files found in the input directory."}), 400
            
            if read_type == "paired":
                fastq_pairs = pair_fastq_files(fastq_files)
                if not fastq_pairs:
                    return jsonify({"error": "No FASTQ files found in the input directory."}), 400
                unpaired = [pair for pair in fastq_pairs if pair[1] is None]
                if unpaired:
                    # Warning will be shown in status updates
                    pass
            else:
                fastq_pairs = [(f, None) for f in fastq_files]

        # Generate job ID
        with job_lock:
            global job_counter
            job_counter += 1
            job_id = str(job_counter)
        
        # Start processing in background thread
        thread = threading.Thread(
            target=process_files_async,
            args=(job_id, genome_dir, output_dir, threads, fastq_pairs)
        )
        thread.daemon = True
        thread.start()
        
        return jsonify({"job_id": job_id, "total_files": len(fastq_pairs)})

    return render_template("index.html")

@app.route("/status/<job_id>")
def get_status(job_id):
    """SSE endpoint for real-time status updates"""
    def generate():
        last_message_count = 0
        while True:
            if job_id in processing_status:
                status = processing_status[job_id]
                current_count = len(status['messages'])
                
                # Send new messages
                if current_count > last_message_count:
                    for i in range(last_message_count, current_count):
                        yield f"data: {json.dumps({'type': 'message', 'text': status['messages'][i]})}\n\n"
                    last_message_count = current_count
                
                # Send status update
                yield f"data: {json.dumps({'type': 'status', 'status': status['status'], 'current': status['current_file'], 'total': status['total_files'], 'multiqc_html': status.get('multiqc_html')})}\n\n"
                
                if status['status'] in ['completed', 'error']:
                    break
            else:
                yield f"data: {json.dumps({'type': 'error', 'text': 'Job not found'})}\n\n"
                break
            
            time.sleep(0.5)
    
    return Response(generate(), mimetype='text/event-stream')

@app.route("/terminal/<job_id>")
def get_terminal(job_id):
    """Get terminal output for a job"""
    if job_id in terminal_output:
        return jsonify({"output": terminal_output[job_id]})
    return jsonify({"output": []})

@app.route("/multiqc")
def serve_multiqc():
    """Serve MultiQC HTML file"""
    filepath = request.args.get('file')
    if not filepath or not os.path.exists(filepath):
        return "File not found", 404
    
    directory = os.path.dirname(filepath)
    filename = os.path.basename(filepath)
    return send_from_directory(directory, filename)

if __name__ == "__main__":
    app.run(debug=True, threaded=True)
