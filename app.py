from flask import Flask, render_template, request, redirect, url_for, flash
import os
import subprocess
from werkzeug.utils import secure_filename

app = Flask(__name__)
app.secret_key = 'supersecretkey'

# Fixed STAR and MultiQC executable paths
STAR_EXECUTABLE = "/home/compbio001/anaconda3/envs/STAR/bin/STAR-avx2"
MULTIQC_EXECUTABLE = "/home/compbio001/anaconda3/envs/multiqc/bin/multiqc"

UPLOAD_FOLDER = 'uploads'
os.makedirs(UPLOAD_FOLDER, exist_ok=True)

@app.route("/", methods=["GET", "POST"])
def index():
    if request.method == "POST":
        genome_dir = request.form.get("genome_dir").strip()
        output_dir = request.form.get("output_dir").strip()
        threads = request.form.get("threads", "8").strip()
        mode = request.form.get("mode")

        if not os.path.exists(genome_dir):
            flash("❌ Invalid genome directory path!", "error")
            return redirect(url_for("index"))

        os.makedirs(output_dir, exist_ok=True)

        fastq_files = []

        if mode == "single":
            uploaded_file = request.files.get("fastq_file")
            if uploaded_file and uploaded_file.filename.endswith(".fastq"):
                filename = secure_filename(uploaded_file.filename)
                file_path = os.path.join(UPLOAD_FOLDER, filename)
                uploaded_file.save(file_path)
                fastq_files = [file_path]
            else:
                flash("❌ Please upload a valid FASTQ file.", "error")
                return redirect(url_for("index"))

        elif mode == "batch":
            input_dir = request.form.get("input_dir").strip()
            if not os.path.exists(input_dir):
                flash("❌ Invalid input directory path!", "error")
                return redirect(url_for("index"))
            fastq_files = [os.path.join(input_dir, f) for f in os.listdir(input_dir) if f.endswith(".fastq")]

        for input_path in fastq_files:
            filename = os.path.basename(input_path)
            sample_name = os.path.splitext(filename)[0]
            output_prefix = os.path.join(output_dir, sample_name)

            star_command = [
                STAR_EXECUTABLE,
                "--runThreadN", threads,
                "--genomeDir", genome_dir,
                "--readFilesIn", input_path,
                "--outFileNamePrefix", output_prefix,
                "--outSAMtype", "BAM", "SortedByCoordinate",
                "--quantMode", "GeneCounts"
            ]

            try:
                subprocess.run(star_command, check=True)
                flash(f"✅ Processed {filename} successfully.", "success")
            except subprocess.CalledProcessError as e:
                flash(f"❌ Error processing {filename}: {str(e)}", "error")

        # Run MultiQC on the output directory and save results in the same directory
        try:
            subprocess.run([MULTIQC_EXECUTABLE, output_dir, "-o", output_dir], check=True)
            flash("📊 MultiQC analysis completed. Please check the QC reports in the output directory.", "success")
        except subprocess.CalledProcessError as e:
            flash(f"❌ Error running MultiQC: {str(e)}", "error")

        return redirect(url_for("index"))

    return render_template("index.html")

if __name__ == "__main__":
    app.run(debug=True, threaded=True)
