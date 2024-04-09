# For the output folder NAD of the MD calculation job of the MD folder, do the data analysis

## Step1：Using convert.py, read the NAD output folder and generate a new SD_basis folder.
          python convert.py (can run in terminal)

## Step2: Generate avg_deco.txt and other output files for different surface hopping methods using decoherence.py
          python decoherence.py (suggest run with sbtach submit.slm)

## Step3: Finally, the tutorial code is used to visualize the data to analyze the excited state decay process. 
