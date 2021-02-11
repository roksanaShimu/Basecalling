# Basecalling
This repository is an implementation of Basecalling for DNA Sequencing. Basecalling is essentially a sequence labeling problem: a time-series is
converted to a sequence of letters from the alphabet A; C; G; T. The aim of this research is
to make the conversion faster and more power efficient to enable embedded operations. I implement a GPU-enhanced sequence analyzer capable of operating in a streaming fashion. 

![block](/basescaller.png)

![viterbi](/viterbi.png)

[Paper](https://ieeexplore.ieee.org/abstract/document/8052869/) 
```
Hossain, R., Mittmann, R., Ghafar-Zadeh, E., Messier, G. G., & Magierowski, S. (2017, August). GPU base calling for DNA strand sequencing. In 2017 IEEE 60th International Midwest Symposium on Circuits and Systems (MWSCAS) (pp. 96-99). IEEE.
```
