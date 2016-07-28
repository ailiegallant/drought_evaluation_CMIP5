

def hello():
	print("hello world!")
	return
	
def numberop(number, op):
	
	if op == "square":
		return square(number)
	elif op == "double":
		return double(number)
	else:
	    print(op+" not found, returning number...")
	    return number

def square(number):
	return number*number
	
def double(number):
	"""This function doubles the input number"""
	return number*2