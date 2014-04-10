#! /usr/bin/python

import re
# Converts the output of selective Mathematica scripts to efficient, syntactically correct C++ code

# input: to_replace.txt, output: replaced.txt

# Set the Mathematica output display mode to Standard Form

text = open('to_replace.txt', 'r').read()


#	Replacement rules
text = re.sub('\n', '', text)
text = re.sub(r'[ ]{2,}', '', text)
text = re.sub(r'(\([a-z]*)\^\\\[Prime\]\\\[Prime\]\)\[t\] \-\>', r'\n\1":\n', text)
text = re.sub(r'Derivative\[1\]\[([a-z]*)\]\[t\]', r'get_d_\1()', text)
text = re.sub(r'([a-z]*)\[t\]', r'get_\1()', text)
text = re.sub(r'\\\[Pi\]', 'M_PI', text)
text = re.sub(r'(([A-Za-z]|[_\(\)])*)\^2', r'square(\1)', text)
text = re.sub(r'(([A-Za-z]|[_\(\)])*)\^3', r'cube(\1)', text)
text = re.sub(r'([^ \t\r\f\v\-\+\\\*\-\>]+)[ \t\r\f\v]+([^ \t\r\f\v\-\+\\\*\-\>]+)', r'\1 * \2', text) # a b c d -> a * b c * d
text = re.sub(r'([^ \t\r\f\v\-\+\\\*\-\>]+)[ \t\r\f\v]+([^ \t\r\f\v\-\+\\\*\-\>]+)', r'\1 * \2', text) # a * b c * d -> a * b * c * d
text = re.sub(r'Sin', 'sin', text)
text = re.sub(r'Cos', 'cos', text)
text = re.sub(r'([0-9]+)', r'\1.0', text)
text = re.sub(r'\[', '(', text)
text = re.sub(r'\]', ')', text)

open('replaced.txt', 'w').write(text)
