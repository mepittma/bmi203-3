# example

[![Build
Status](https://travis-ci.org/mepittma/bmi203-3.svg?branch=master)](https://travis-ci.org/mepittma/bmi203-3)

Example python project with testing.

## usage

To use the package, first make a new conda environment and activate it

```
conda create -n align_env python=3
source activate align_env
```

then run

```
conda install --yes --file requirements.txt
```

to install all the dependencies in `requirements.txt`. Then the package's
main function (located in `align/__main__.py`) can be run as follows

```
python -m align
```

## testing

Testing is as simple as running

```
python -m pytest
```

from the root directory of this project.
