# uavrt_noise
A repo for noise assessment on our systems.

To install the necessary pacakges, navigate to this repo in terminal then:
```
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```
Once installed and for subsequent runs, you need to source the virtual environment before you run the progam with:
```
source .venv/bin/activate
```

To run:
```
python 3 noisesweeper.py "DescriptiveNameofNoiseTest"
```

To process results, be sure to fist modify the lines in `noiseresultprocessor.py` below the line that says `# MODIFY THE FOLLOWING LINES OF CODE FOR YOUR SPECIFIC TEST`. Then run
```
python3 noiseresultprocessor.py
```

This will produce two plots, one is the raw results power spectral density, the other is a 3750,000 Hz moving average result of the power spectral density. 

