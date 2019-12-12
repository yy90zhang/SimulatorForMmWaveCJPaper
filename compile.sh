#!/bin/sh

mpic++ -o simulator.out *.cpp -std=c++11

mv simulator.out ../..
