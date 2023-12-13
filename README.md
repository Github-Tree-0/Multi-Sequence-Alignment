# Multi-Sequence-Alignmen

### Authors:
ClustalW (Rakshitha Kiran - RakKiran09)\
Feng and Doolittle (Jipeng Lyu - Github-Tree-0)

### Environment
To run our code, please prepare the environment first by running
```
pip install matplotlib numpy pandas scipy
```

### ClustalW
This folder contains all the implementations of the ClustalW algorithm. To run the code, please run
```
python ./main.py
```
in the `./ClustalW/` folder. You can change the input sequences in `./ClustalW/main.py`, 
but we recommand you refer to our test code `./test.ipynb` to play with it.

### Feng and Doolittle
The folder `./Feng_Doolittle` contains all the implementations of the Feng and Doolittle algorithm. 
In `./Feng_Doolittle/algorithm.py`, I defined a MSA class, and refered to the class in `./Feng_Doolittle/test_space.ipynb`.
To play with the code, please use `./Feng_Doolittle/test_space.ipynb`, or you can refer to `./test.ipynb`

### Experiment
We also provide our experiment code on both synthetic and real dataset in `./test.ipynb`. The synthetic dataset generation function is also there.
