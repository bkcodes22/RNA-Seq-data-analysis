import subprocess
import inquirer


ascii_art = r"""
  __          __  _                            _          _____  ______ _____ ______                
 \ \        / / | |                          | |        |  __ \|  ____/ ____|  ____|               
  \ \  /\  / ___| | ___ ___  _ __ ___   ___  | |_ ___   | |  | | |__ | |  __| |__   __ _ ___ _   _ 
   \ \/  \/ / _ | |/ __/ _ \| '_ ` _ \ / _ \ | __/ _ \  | |  | |  __|| | |_ |  __| / _` / __| | | |
    \  /\  |  __| | (_| (_) | | | | | |  __/ | || (_) | | |__| | |___| |__| | |___| (_| \__ | |_| |
     \/  \/ \___|_|\___\___/|_| |_| |_|\___|  \__\___/  |_____/|______\_____|______\__,_|___/\__, |
                                                                                              __/ |
                                                                                             |___/ 
                                                                                                                                                                                                                                                                                                                                                                                                                 
"""
print(ascii_art)

def run_fastq_processing():
    """Run the Python app.py for Fastq Processing"""
    print("Running Fastq Processing...")
    ascii_art = r"""
                                .-----.
                               / .===. \
                               \/ 6 6 \/
                               ( \___/ )
  _________________________ooo__\_____/______________________________
 /                                                                   \
| Let's get the counts file first, copy the below IP to your browser! |
 \_______________________________________ooo_________________________/
                                |  |  |
                                |_ | _|
                                |  |  |
                                |__|__|
                                /-'Y'-\
                               (__/ \__)                                                                                                                                                                                                                                                                                                                                                                                               
    """
    print(ascii_art)
    subprocess.run(["python", "app.py"])

def run_deg_analysis():
    """Run the R app.R for DEG analysis"""
    print("Running DEG Analysis...")
    ascii_art = r"""

               /\             /\
              |`\\_,--="=--,_//`|
              \ ."  :'. .':  ". /
             ==)  _ :  '  : _  (==
               |>/O\   _   /O\<|
               | \-"~` _ `~"-/ |
              >|`===. \_/ .===`|<
        .-"-.   \==='  |  '===/   .-"-.
.------{'. '`}---\,  .-'-.  ,/---{.'. '}-----.
 )     `"---"`     `~-===-~`     `"---"`    (
(  Streamlining RNA-Seq, One DEG at a Time!  )
 )                                          (
'--------------------------------------------'
"""
    print(ascii_art)
    subprocess.run(["Rscript", "shiny_app/app.R"])

def main():
    # Create a list of questions for the interactive CLI
    questions = [
        inquirer.List(
            'task',
            message="What would you like to do?",
            choices=['Fastq Processing', 'DEG Analysis'],
        ),
    ]
    
    # Get the answer from the user
    answers = inquirer.prompt(questions)
    
    # Check the user's selection and call the corresponding function
    if answers['task'] == 'Fastq Processing':
        run_fastq_processing()
    elif answers['task'] == 'DEG Analysis':
        run_deg_analysis()

if __name__ == "__main__":
    main()
