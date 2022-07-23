# Mathematical formulations for graph burning
This repository contains two mathematical formulations for the graph burning problem:

| Name  | Description | Binarry Variables | Integer variables | Constraints |
| ------------- | ------------- | ------------- | ------------- | ------------- |
| ILP  | integer linear program  | O(nU) | 0 | O(\|E\|U) |
| CSP1  | constraint satisfaction problem 1 | O(nB) | 0 | O(\|E\|U) |
| CSP2  | constraint satisfaction problem 2 | O(n^2) | O(n) | O(n^2) |

CSP1+BS and CSP2+BS consist in adding each formulation within a bnary search. This way, each formulation is executed log n times with different guesses B on b(G).

To execute the implemented formulations you need to install Gurobi.

## Install gurobipy:

Source: https://www.gurobi.com/gurobi-and-anaconda-for-windows/

### Step one: Download and install Anaconda

Gurobi supports Python 2.7 and 3.7 for Windows. However, to run our code install Python 3.X. Please choose the version of Anaconda you wish to download (the download will start automatically):

Once the download is complete, click on it to run the installer.

### Step two: Install Gurobi into Anaconda

The next step is to install the Gurobi package into Anaconda. You do this by first adding the Gurobi channel into your Anaconda platform and then installing the gurobi package from this channel.

From an Anaconda terminal issue the following command to add the Gurobi channel to your default search list:

```
$ conda config --add channels http://conda.anaconda.org/gurobi
```

Now issue the following command to install the Gurobi package:

```
$ conda install gurobi
```

You can remove the Gurobi package at any time by issuing the command:

```
$ conda remove gurobi
```

### Step three: Install a Gurobi License

The third step is to install a Gurobi license (if you haven’t already done so).

You are now ready to use Gurobi from within Anaconda. Your next step is to launch either the Spyder IDE or Jupyter Notebook.


## Running the implemented formulations

Change "folder" and "folder_dataset" in the code to point to your dataset. Then, execute the code in the Anaconda prompt.

```
$ python ILP.py
```
or
```
$ python CSP1+BS.py
```
or
```
$ python CSP2+BS.py
```

The returned optimal burning sequence have vertices enumerated from 0 to n.

## Optionally, you can verify and draw your solutions

Use the Jupyter notebook 'Verify and Draw' to get images like the following:

![image](https://user-images.githubusercontent.com/9120537/166972154-37f314df-d90e-4152-8443-9d64fe83dd52.png)

![image](https://user-images.githubusercontent.com/9120537/166973826-10c61d4a-d922-4369-b7ac-c2c902f9ba51.png)

![image](https://user-images.githubusercontent.com/9120537/166975421-81a684a5-1940-4839-bb00-87ef163099fb.png)

## Anyway, What is graph burning?

[![Click me](https://img.youtube.com/vi/oxX8ONv_8Pk/0.jpg)](https://www.youtube.com/watch?v=oxX8ONv_8Pk)
