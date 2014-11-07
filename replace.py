#! /usr/bin/python

import re
import string
# Converts the output of selective Mathematica scripts to efficient, syntactically correct C++ code

# input: to_replace.txt, output: replaced.txt

# Set the Mathematica output display mode to Standard Form

text = open('LagrangianSolutions.txt', 'r').read()
out = open('dynein_motion_functions.cpp', 'w')

out.write("#include \"dynein_struct.h\"\n\n")

#	Replacement rules
text = re.sub(r'Derivative\(bla\(t\), t\)', 'get_d_bra()', text)
text = re.sub(r'Derivative\(mla\(t\), t\)', 'get_d_bra()', text)
text = re.sub(r'Derivative\(mra\(t\), t\)', 'get_d_bra()', text)
text = re.sub(r'Derivative\(bra\(t\), t\)', 'get_d_bra()', text)

text = re.sub(r'bla\(t\)', 'get_bla()', text)
text = re.sub(r'mla\(t\)', 'get_mla()', text)
text = re.sub(r'mra\(t\)', 'get_mra()', text)
text = re.sub(r'bra\(t\)', 'get_bra()', text)

text = re.sub(r'[ ]+', '', text)

#text = re.sub(r'([a-zA-Z]+)\^2', r'square(\1)', text)
#text = re.sub(r'([a-zA-Z]+)\^3', r'cube(\1)', text)

idx = 0
while string.find(text, ")**2") != -1: # Because nested paren regex (\((?:[^\(\)]+|(?R))*\))(?:\^2)? doesn't work for this pattern :( 
	idx = string.find(text, ")**2")
	i = idx
	v = 1
	while (v != 0):
		i = i - 1			# Iterate backwards until matching parenthesis is found
		if text[i] == ")":
			v = v + 1
		elif text[i] == "(":
			v = v - 1
	if text[i-3:i] == "cos":
		text = text[:i-3] + "square(" + text[i-3:idx+1] + ")" + text[idx+4:]
			
	elif text[i-3:i] == "sin":
		text = text[:i-3] + "square(" + text[i-3:idx+1] + ")" + text[idx+4:]
		
	elif text[i-9:i-5] == "get_":
		text = text[:i-9] + "square(" + text[i-9:idx+1] + ")" + text[idx+4:]
		
	else:
		text = text[:i] + "square" + text[i:idx+1] + text[idx+4:]
	
idx = 0
while string.find(text, ")**3") != -1:
	idx = string.find(text, ")**3")
	i = idx
	v = 1
	while (v != 0):
		i = i - 1
		if text[i] == ")":
			v = v + 1
		elif text[i] == "(":
			v = v - 1
	if text[i-3:i] == "cos":
		text = text[:i-3] + "cube(" + text[i-3:idx+1] + ")" + text[idx+4:]
		
	elif text[i-3:i] == "sin":
		text = text[:i-3] + "cube(" + text[i-3:idx+1] + ")" + text[idx+4:]
		
	elif text[i-9:i-5] == "get_":
		text = text[:i-9] + "cube(" + text[i-9:idx+1] + ")" + text[idx+4:]
		
	text = text[:i] + "cube" + text[i:idx+1] + text[idx+4:]
	
text = re.sub(r'([A-Za-z]*)\*\*2', r'square(\1)', text)
text = re.sub(r'([A-Za-z]*)\*\*3', r'cube(\1)', text)
text = re.sub(r'([^ \t\r\f\v\-\+\\\*\-\>]+)[ \t\r\f\v]+([^ \t\r\f\v\-\+\\\*\-\>]+)', r'\1 * \2', text) # a b c d -> a * b c * d
text = re.sub(r'([^ \t\r\f\v\-\+\\\*\-\>]+)[ \t\r\f\v]+([^ \t\r\f\v\-\+\\\*\-\>]+)', r'\1 * \2', text) # a * b c * d -> a * b * c * d
#text = re.sub(r'([0-9]+)', r'\1.0', text)

text = re.sub(r'\{Derivative\(get_bla\(\),t,t\):([^\n]*),\n', r'double Dynein::get_dd_bla() {\nreturn \1;\n}\n\n', text)
text = re.sub(r'Derivative\(get_mla\(\),t,t\):([^\n]*),\n', r'double Dynein::get_dd_mla() {\nreturn \1;\n}\n\n', text)
text = re.sub(r'Derivative\(get_mra\(\),t,t\):([^\n\}]*)\}\n', r'double Dynein::get_dd_mra() {\nreturn \1;\n}\n\n', text)
text = re.sub(r'Derivative\(get_bra\(\),t,t\):([^\n]*),\n', r'double Dynein::get_dd_bra() {\nreturn \1;\n}\n\n', text)

out.write(text)
