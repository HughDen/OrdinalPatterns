This project contains code for a recent [paper](https://arxiv.org/abs/1907.07172) I wrote about ordinal patterns that occur in random walks. All of the notebooks are Python 3. With the exception of a few print statements, I expect most of this code works in Python 2. 

The notebooks contained in the Calculation_notebooks folders perform some calculation related to the paper.


[Alcoves_in_the_polytope.ipynb](https://github.com/HughDen/OrdinalPatterns/blob/master/Calculation_notebooks/Alcoves_in_the_polytope.ipynb) uses a brute force breadth-first search algorithm to generate all addresses of alcoves for the polytope defined in the paper.

[Uniform_dictionary.ipynb](https://github.com/HughDen/OrdinalPatterns/blob/master/Calculation_notebooks/Uniform_dictionary.ipynb) uses a brute force method that counts alcoves to compute ordinal pattern probabilities for random walks whose steps are uniform.

[Almost_consecutive_totals.ipynb](https://github.com/HughDen/OrdinalPatterns/blob/master/Calculation_notebooks/Almost_consecutive_totals.ipynb) calculates the probability that a random walk with symmetric steps produces an almost consecutive permutation as an ordinal pattern.

The notebooks with in the Test_notebooks folder run simulations and compare to the exact results we expect using the theorems of the paper.

Right now, the project just contains code that a reader of the paper might expect. It is the code I used to explore some things during the research process. There is currently no class hierarchy. It's all just functions. Any suggestions or comments (Hugh.Denoncourt@Colorado.edu) are welcome.