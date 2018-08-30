#!/bin/bash

ifort -c stieltjes.f90 imaging.f90 interp.f90 pythag.f sort2.f90 tql2.f
ifort stieltjes.o imaging.o interp.o pythag.o sort2.o tql2.o -o stieltjes
