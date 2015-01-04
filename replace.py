#! /usr/bin/python

import re
import string
# Converts the output of selective Mathematica scripts to efficient, syntactically correct C++ code

# input: DyneinBrownianBothboundSolutionsUnsimplified.txt, DyneinBrownianLeftboundSolutionsUnsimplified.txt, DyneinBrownianRightboundSolutionsUnsimplified.txt, output: replaced.txt

# Set the Mathematica output display mode to Standard Form

text_both =  open('DyneinBrownianBothboundSolutionsUnsimplified.txt', 'r').read()
text_left =  open('DyneinBrownianLeftboundSolutionsUnsimplified.txt', 'r').read()
text_right = open('DyneinBrownianRightboundSolutionsUnsimplified.txt', 'r').read()

text = text_left + text_right + text_both
out = open('dynein_motion_functions.cpp', 'w')

out.write("#include \"dynein_struct.h\"\n\n")

#	Replacement rules
text = re.sub('\n', '', text)															# Remove all line breaks
text = re.sub(r'[ ]{2,}', '', text)			 											# Delete all large spaces
text = re.sub(r'(\([a-z]*)\^\\\[Prime\]\\\[Prime\]\)\[t\] \-\>', r'\n\1":\n', text)		# Change 'Derivative[1]...' to 'get_d_...'
text = re.sub(r'Derivative\[1\]\[([a-z]*)\]\[t\]', r'get_d_\1()', text)					
text = re.sub(r'([a-z]*)\[t\]', r'get_\1()', text)										# Change 'bla[t]' etc to 'get_bla(t)'
text = re.sub(r'\\\[Pi\]', 'M_PI', text)												# Write pi in C
text = re.sub(r'Sin', 'sin', text)														
text = re.sub(r'Cos', 'cos', text)
text = re.sub(r'\[', '(', text)															# Change Mathematica brackets to C parens
text = re.sub(r'\]', ')', text)
text = re.sub(r'([a-zA-Z]+)\^2', r'square(\1)', text)									# Change ^2/^3 to local square/cube functions for single-variable expressions
text = re.sub(r'([a-zA-Z]+)\^3', r'cube(\1)', text)

idx = 0
while string.find(text, ")^2") != -1: # Change ^2/^3 to local square/cube functions for multi-variable expressions	
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

text = re.sub(r'\{\{Derivative\(2\.0\)\(bla\)get_\(\) - \>([^,\n]*),', r'bla: \1\n', text)		# Build temp structure to hold different expressions
text = re.sub(r'Derivative\(2\.0\)\(mla\)get_\(\) - \>([^,\n]*),', r'mla: \1\n', text)
text = re.sub(r'Derivative\(2\.0\)\(mra\)get_\(\) - \>([^,\n]*),', r'mra: \1\n', text)
text = re.sub(r'Derivative\(2\.0\)\(bra\)get_\(\) - \>([^\}\n]*)\}\}', r'bra: \1\n', text)

text = re.sub(r' \* ', '', text)

text = re.sub(r'bla: (.*)\nmla: (.*)\nmra: (.*)\nbra: (.*)\nbla: (.*)\nmla: (.*)\nmra: (.*)\nbra: (.*)\nbla: (.*)\nmla: (.*)\nmra: (.*)\nbra: (.*)\n',		# Create acceleration functions
r'bla: \1\nmla: \2\nmra: \3\nbra: \4\nbla: \5\nmla: \6\nmra: \7\nbra: \8\nbla: \9\nmla: \10\nmra: \11\nbra: \12\n' + 
r'double Dynein::get_dd_bla() {\n\tif (state == BOTHBOUND)\n\t\treturn \1;\n\telse if (state == LEFTBOUND)\n\t\treturn \5;\n\telse\n\t\treturn \9;\n}\n\n', text)

text = re.sub(r'bla: (.*)\nmla: (.*)\nmra: (.*)\nbra: (.*)\nbla: (.*)\nmla: (.*)\nmra: (.*)\nbra: (.*)\nbla: (.*)\nmla: (.*)\nmra: (.*)\nbra: (.*)\n',
r'bla: \1\nmla: \2\nmra: \3\nbra: \4\nbla: \5\nmla: \6\nmra: \7\nbra: \8\nbla: \9\nmla: \10\nmra: \11\nbra: \12\n' + 
r'double Dynein::get_dd_mla() {\n\tif (state == BOTHBOUND)\n\t\treturn \2;\n\telse if (state == LEFTBOUND)\n\t\treturn \6;\n\telse\n\t\treturn \10;\n}\n\n', text)

text = re.sub(r'bla: (.*)\nmla: (.*)\nmra: (.*)\nbra: (.*)\nbla: (.*)\nmla: (.*)\nmra: (.*)\nbra: (.*)\nbla: (.*)\nmla: (.*)\nmra: (.*)\nbra: (.*)\n',
r'bla: \1\nmla: \2\nmra: \3\nbra: \4\nbla: \5\nmla: \6\nmra: \7\nbra: \8\nbla: \9\nmla: \10\nmra: \11\nbra: \12\n' + 
r'double Dynein::get_dd_mra() {\n\tif (state == BOTHBOUND)\n\t\treturn \3;\n\telse if (state == LEFTBOUND)\n\t\treturn \7;\n\telse\n\t\treturn \11;\n}\n\n', text)

text = re.sub(r'bla: (.*)\nmla: (.*)\nmra: (.*)\nbra: (.*)\nbla: (.*)\nmla: (.*)\nmra: (.*)\nbra: (.*)\nbla: (.*)\nmla: (.*)\nmra: (.*)\nbra: (.*)\n',
r'bla: \1\nmla: \2\nmra: \3\nbra: \4\nbla: \5\nmla: \6\nmra: \7\nbra: \8\nbla: \9\nmla: \10\nmra: \11\nbra: \12\n' + 
r'double Dynein::get_dd_bra() {\n\tif (state == BOTHBOUND)\n\t\treturn \4;\n\telse if (state == LEFTBOUND)\n\t\treturn \8;\n\telse\n\t\treturn \12;\n}\n\n', text)

text = re.sub(r'bla: (.*)\nmla: (.*)\nmra: (.*)\nbra: (.*)\nbla: (.*)\nmla: (.*)\nmra: (.*)\nbra: (.*)\nbla: (.*)\nmla: (.*)\nmra: (.*)\nbra: (.*)\n', '', text)		# Remove temp structure

out.write(text)
