#######################################
# ambre.mk
# Default options for ambre computer
#######################################
CC=gcc
LIBSLOCAL=-L/opt/homebrew/opt/lapack/lib -llapack -L/opt/homebrew/opt/openblas/lib -lopenblas -L/usr/lib -lm
INCLUDEBLASLOCAL=-I/usr/include -I/opt/homebrew/opt/lapack/include -I/opt/homebrew/opt/openblas/include
OPTCLOCAL=-fPIC -march=native
