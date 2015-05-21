#! /usr/bin/python

import re
import string
# Converts the output of Mathematica scripts (in Standard Form) in Mathematica/ to C++ code

# input: LeftboundSolutions.m
# output: dynein_motion_functions.cpp

text =  open('Motion_Equations/Leftbound_Solution.txt', 'r').read()

out = open('dynein_motion_functions.cpp', 'w')

out.write("#include \"dynein_struct.h\"\n\n")

#line = 0
#while (string.find(text, "\n")) != -1:												# Annotate output .cpp with line numbers of original file
	#idx = string.find(text, "\n")
	#text = text[:idx] + "/* " + str(int(line)) + " */" + text[idx+2:]
	#line += 1

#	Replacement rules

#text = re.sub(r'{{([^{}]*)},{([^{}]*)},\n{([^{}]*)},\n{([^{}]*)},\n{([^{}]*)},\n{([^{}]*)},\n{([^{}]*)},\n{([^{}]*)}}', r'TVBL:\n\1\n\n TVBL:\n\2\n\n TVBL:\n\3\n\n TVBL:\n\4\n\n TVBL:\n\5\n\n TVBL:\n\6\n\n TVBL:\n\7\n\n TVBL:\n\8\n\n ', text)
text = re.sub(r'\s', '', text)
text = re.sub(r'{{([^{}]*)},{([^{}]*)},{([^{}]*)},{([^{}]*)},{([^{}]*)},{([^{}]*)},{([^{}]*)},{([^{}]*)}}', r'TVBL: \1\nTVML: \2\nTVMR: \3\nTVBR: \4\nLBL: \5\nLML: \6\nLMR: \7\nLBR: \8\n', text)

print "replace.py: Changing Mathematica symbols to C symbols"

text = re.sub('Pi', 'M_PI', text)														# Write pi in C
text = re.sub(r'Sin', 'sin', text)														
text = re.sub(r'Cos', 'cos', text)
text = re.sub(r'\[', '(', text)															# Change Mathematica brackets to C parens
text = re.sub(r'\]', ')', text)
text = re.sub(r'([a-zA-Z]+)\^2', r'square(\1)', text)									# Change ^2/^3 to local square/cube functions for single-variable expressions
text = re.sub(r'([a-zA-Z]+)\^3', r'cube(\1)', text)
text = re.sub(r'([a-zA-Z]+)\^4', r'fourth(\1)', text)
text = re.sub(r'([a-zA-Z]+)\^5', r'fifth(\1)', text)
text = re.sub(r'g', r'1/g', text)

text = re.sub(r'Dxbl', r' -ls*sin(get_bla())', text)
text = re.sub(r'Dxml', r' -lt*sin(get_mla())', text)
text = re.sub(r'Dxmr', r' -lt*sin(get_mra()-M_PI)', text)
text = re.sub(r'Dxbr', r' -ls*sin(get_bra()-M_PI)', text)

text = re.sub(r'Dybl', r' -ls*cos(get_bla())', text)
text = re.sub(r'Dyml', r' -lt*cos(get_mla())', text)
text = re.sub(r'Dymr', r' -lt*cos(get_mra()-M_PI)', text)
text = re.sub(r'Dybr', r' -ls*cos(get_bra()-M_PI)', text)

text = re.sub(r'Tbl', r'get_bla()', text)
text = re.sub(r'Tml', r'get_mla()', text)
text = re.sub(r'Tmr', r'get_mra()', text)
text = re.sub(r'Tbr', r'get_bra()', text)

print "replace.py: Changing Mathematica variable names to C variable names"
text = re.sub('Lt', 'lt', text)															# Convert from Mathematica variable naming scheme to C variable naming scheme
text = re.sub('Ls', 'ls', text)
text = re.sub(r'F([mbrl]{2})x', r'f\1x', text)
text = re.sub(r'F([mbrl]{2})y', r'f\1y', text)

print "replace.py: Changing Mathematica squares / cubes to C squares / cubes"
idx = 0
while string.find(text, ")^2") != -1:													 # Change ^2/^3 to local square/cube functions for multi-variable expressions	
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
		continue
			
	elif text[i-3:i] == "sin":
		text = text[:i-3] + "square(" + text[i-3:idx+1] + ")" + text[idx+3:]
		continue
		
	elif text[i-9:i-5] == "get_":
		text = text[:i-9] + "square(" + text[i-9:idx+1] + ")" + text[idx+3:]
		continue
		
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
		continue
		
	elif text[i-3:i] == "sin":
		text = text[:i-3] + "cube(" + text[i-3:idx+1] + ")" + text[idx+3:]
		continue
		
	elif text[i-9:i-5] == "get_":
		text = text[:i-9] + "cube(" + text[i-9:idx+1] + ")" + text[idx+3:]
		continue
		
	else:
		text = text[:i] + "cube" + text[i:idx+1] + text[idx+3:]

"""

print "replace.py: Adding multiplication symbols"
text = re.sub(r'([^ \t\r\f\v\-\+\\\*\-\>]+)[ \t\r\f\v]+([^ \t\r\f\v\-\+\\\*\-\>]+)', r'\1 * \2', text) # a b c d -> a * b c * d
text = re.sub(r'([^ \t\r\f\v\-\+\\\*\-\>]+)[ \t\r\f\v]+([^ \t\r\f\v\-\+\\\*\-\>]+)', r'\1 * \2', text) # a * b c * d -> a * b * c * d
text = re.sub(r'([0-9]+)', r'\1.0', text)
text = re.sub(r'([\+-]{1})', r'\1 ', text)

print "replace.py: Creating CPP file"
text = re.sub(r'BLA_DERIVATIVE_2\.0:([^,\n]*),', r'bla: \1\n', text)								   # Build temp structure to hold different expressions
text = re.sub(r'MLA_DERIVATIVE_2\.0:([^,\n]*),', r'mla: \1\n', text)
text = re.sub(r'MRA_DERIVATIVE_2\.0:([^,\n]*),', r'mra: \1\n', text)
text = re.sub(r'BRA_DERIVATIVE_2\.0:([^\}\n]*)\}\}', r'bra: \1\n', text)

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
"""

text = re.sub(r'TVBL: (.*)\nTVML: (.*)\nTVMR: (.*)\nTVBR: (.*)\nLBL: (.*)\nLML: (.*)\nLMR: (.*)\nLBR: (.*)\n',
	r'double Dynein::get_d_bla() {\n\tif (state == LEFTBOUND)\n\t\treturn \1;\n\telse if (state == RIGHTBOUND)\n\t\treturn 0;\n\telse\n\t\treturn 0;\n}\n\n' + 
	r'double Dynein::get_d_mla() {\n\tif (state == LEFTBOUND)\n\t\treturn \2;\n\telse if (state == RIGHTBOUND)\n\t\treturn 0;\n\telse\n\t\treturn 0;\n}\n\n' + 
	r'double Dynein::get_d_mra() {\n\tif (state == LEFTBOUND)\n\t\treturn \3;\n\telse if (state == RIGHTBOUND)\n\t\treturn 0;\n\telse\n\t\treturn 0;\n}\n\n' + 
	r'double Dynein::get_d_bra() {\n\tif (state == LEFTBOUND)\n\t\treturn \4;\n\telse if (state == RIGHTBOUND)\n\t\treturn 0;\n\telse\n\t\treturn 0;\n}\n\n' + 
	r'double Dynein::get_ls_force() {\n\tif (state == LEFTBOUND)\n\t\treturn \5;\n\telse if (state == RIGHTBOUND)\n\t\treturn 0;\n\telse\n\t\treturn 0;\n}\n\n' + 
	r'double Dynein::get_lt_force() {\n\tif (state == LEFTBOUND)\n\t\treturn \6;\n\telse if (state == RIGHTBOUND)\n\t\treturn 0;\n\telse\n\t\treturn 0;\n}\n\n' + 
	r'double Dynein::get_rt_force() {\n\tif (state == LEFTBOUND)\n\t\treturn \7;\n\telse if (state == RIGHTBOUND)\n\t\treturn 0;\n\telse\n\t\treturn 0;\n}\n\n' + 
	r'double Dynein::get_rs_force() {\n\tif (state == LEFTBOUND)\n\t\treturn \8;\n\telse if (state == RIGHTBOUND)\n\t\treturn 0;\n\telse\n\t\treturn 0;\n}\n\n', text)
	
text = re.sub(r'Fybl', r'get_f_bly()', text)
text = re.sub(r'Fyml', r'get_f_mly()', text)
text = re.sub(r'Fyt',  r'get_f_ty()', text)
text = re.sub(r'Fymr', r'get_f_mry()', text)
text = re.sub(r'Fybr', r'get_f_bry()', text)

text = re.sub(r'Fxbl', r'get_f_blx()', text)
text = re.sub(r'Fxml', r'get_f_mlx()', text)
text = re.sub(r'Fxt',  r'get_f_tx()', text)
text = re.sub(r'Fxmr', r'get_f_mrx()', text)
text = re.sub(r'Fxbr', r'get_f_brx()', text)


text = re.sub(r'Rybl', r'get_r_bly()', text)
text = re.sub(r'Ryml', r'get_r_mly()', text)
text = re.sub(r'Ryt',  r'get_r_ty()', text)
text = re.sub(r'Rymr', r'get_r_mry()', text)
text = re.sub(r'Rybr', r'get_r_bry()', text)

text = re.sub(r'Rxbl', r'get_r_blx()', text)
text = re.sub(r'Rxml', r'get_r_mlx()', text)
text = re.sub(r'Rxt',  r'get_r_tx()', text)
text = re.sub(r'Rxmr', r'get_r_mrx()', text)
text = re.sub(r'Rxbr', r'get_r_brx()', text)

out.write(text)
