# Copyright Oscar Dowson, 2017

dominates(::Maximisation, x, y) = x < y
dominates(::Minimisation, x, y) = x > y
