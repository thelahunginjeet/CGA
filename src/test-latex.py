# simple test script to see if we can get around LaTeX issues and string formatting


from CGAGenerator import CGAGenerator


texHeader = r"""
\documentclass[11pt]{amsart}
\usepackage{geometry}                
\geometry{letterpaper}                   
\usepackage{graphicx}
\usepackage{amssymb}
\usepackage{amsmath}
\usepackage{epstopdf}

\title{Genetic Programming Function List}
\author{CBKB}
\begin{document}
\maketitle
"""

equation = r"""

\begin{equation}
%s
\end{equation}
"""

texFooter = r"""
\end{document}
"""

def main():
	# generate a bunch of trees 
	trees = map(CGAGenerator.generate, range(10))
	latex = [equation%(tree.getLatex().replace('$','\\')) for tree in trees]
	string = [equation%(tree.getString()) for tree in trees]
	fo = open('/Users/brown/Work/Temp/junk.tex','w')
	fo.write(texHeader)
	for i in range(len(latex)):
		fo.write(latex[i])
		fo.write(string[i])
	fo.write(texFooter)
	fo.close()

if __name__ == '__main__':
	main()
