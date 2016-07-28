#THIS MODULE CONTAINS SCRIPTS FOR READING ASCII (TEXT) FILES
#EACH SUBROUTINE IS SEPARATED BY A LINE --------
#-----------------------------------------------------------------------------

#
#Reads a single column ASCII
#
#INPUT
#  file - input file name, including full path, of the ASCII file to be read
#
#

def read_sglcol_ascii(file):
    import numpy as np #import numpy for array manipulation
    
    
    f = open(file, 'r')  #open the input file, FILE
   
    line = f.readlines()  #read individual lines in file
    out = [float(i) for i in line]  #convert lines in list to floating points
    arr = np.array(out)   #convert list to array
    
    f.close()  #close the file

    return arr     #return the array

#-----------------------------------------------------------------------------
#
#Reads a multi-column ASCII
#
#INPUT
#  file - input file name, including full path, of the ASCII file to be read
#
#

