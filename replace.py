#! /usr/bin/python

import re
import string
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
text = re.sub(r'Sin', 'sin', text)
text = re.sub(r'Cos', 'cos', text)
text = re.sub(r'\[', '(', text)
text = re.sub(r'\]', ')', text)
text = re.sub(r'([a-zA-Z]+)\^2', r'square(\1)', text)
text = re.sub(r'([a-zA-Z]+)\^3', r'cube(\1)', text)

idx = 0
while string.find(text, ")^2") != -1: # Because nested paren regex (\((?:[^\(\)]+|(?R))*\))(?:\^2)? doesn't work for this pattern :( 
	idx = string.find(text, ")^2")
	i = idx
	v = 1
	while (v != 0):
		i = i - 1
		if text[i] == ")":
			v = v + 1
		elif text[i] == "(":
			v = v - 1
	if text[i-3:i] == "cos":
		text = text[:i-3] + "square(" + text[i-3:idx+1] + ")" + text[idx+3:]
			
	elif text[i-3:i] == "sin":
		text = text[:i-3] + "square(" + text[i-3:idx+1] + ")" + text[idx+3:]
		
	elif text[i-9:i-5] == "get_":
		text = text[:i-9] + "square(" + text[i-9:idx+1] + ")" + text[idx+3:]
		
	else:
		text = text[:i] + "square" + text[i:idx+1] + text[idx+3:]
	
idx = 0
while string.find(text, ")^3") != -1:
	idx = string.find(text, ")^3")
	i = idx
	v = 1
	while (v != 0):
		i = i - 1
		if text[i] == ")":
			v = v + 1
		elif text[i] == "(":
			v = v - 1
	if text[i-3:i] == "cos":
		text = text[:i-3] + "cube(" + text[i-3:idx+1] + ")" + text[idx+3:]
		
	elif text[i-3:i] == "sin":
		text = text[:i-3] + "cube(" + text[i-3:idx+1] + ")" + text[idx+3:]
		
	elif text[i-9:i-5] == "get_":
		text = text[:i-9] + "cube(" + text[i-9:idx+1] + ")" + text[idx+3:]
		
	text = text[:i] + "cube" + text[i:idx+1] + text[idx+3:]

text = re.sub(r'([^ \t\r\f\v\-\+\\\*\-\>]+)[ \t\r\f\v]+([^ \t\r\f\v\-\+\\\*\-\>]+)', r'\1 * \2', text) # a b c d -> a * b c * d
text = re.sub(r'([^ \t\r\f\v\-\+\\\*\-\>]+)[ \t\r\f\v]+([^ \t\r\f\v\-\+\\\*\-\>]+)', r'\1 * \2', text) # a * b c * d -> a * b * c * d
text = re.sub(r'([0-9]+)', r'\1.0', text)
text = re.sub(r'([\+-]{1})', r'\1 ', text)
text = re.sub(r'(.{130}[^\s]*\s)', r'\1\n\t\t', text)

open('replaced.txt', 'w').write(text)
