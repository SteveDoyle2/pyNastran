# coding: utf-8
"""
This code is designed to help you find what's wrong with an .ipynb
file, so you can hand fix it.  The issue with encoding errors is
that it's always custom.  Finding the characters is the biggest
problem.

An example, the encoding in a file should have been utf-8 and was
actually latin1.  As such, the shape of a matrix (e.g., 100 x 200)
was written with a fancy x character.

I changed the x and it was fine.  After rerunning, the fancy x was
back, but due to the newer version of the Jupyter Notebook (and likely
many bug fixes between when it was made and now), the new file works.
"""
import sys
import codecs
with codecs.open(sys.argv[1], 'r', encoding='latin1') as f:
    lines = f.readlines()

text = ''
for line in lines:
    print(line.rstrip())
    text += line.rstrip() + '\n'

text = text.encode('utf8')
new_file_path = 'new/path'
with open(sys.argv[1] + '.new', 'wb') as f:
    f.write(text)
