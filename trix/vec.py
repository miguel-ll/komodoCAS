from numpy import sqrt,arccos,sign

class Quaternion(object):
	def __init__(self,*args):
		if len(args[0])!=4: self.values = [0,0,0,0]
		else: self.values = [a for a in args[0]]

	def __add__(self,other):
		if isinstance(other, Quaternion):
		    added = [ a + b for a, b in zip(self, other) ]
		elif isinstance(other, (int, float)):
		    added = [a + other for a in self]
		else:
		    raise ValueError("Addition with type {} not supported".format(type(other)))
		return self.__class__(added)
	def __iter__(self):
		return self.values.__iter__()
	def __len__(self):
		return len(self.values)
	def __getitem__(self, key):
		return self.values[key]
	def __repr__(self):
		return str(self.values)

#v = Quaternion([1,2,3,4])
#u = Quaternion([5,6,7,8])
#print(u+v)

class Vector(object):
    def __init__(self, *args):
        if len(args)==0: self.values = [0,0,0]
        else: self.values = [a for a in args[0]]

    def norm(self):
        return sqrt(sum(x*x for x in self))

    def unit(self):
        # Returns a normalized unit vector
        norm = self.norm()
        normed = [x/norm for x in self]
        return self.__class__(normed)

    def __add__(self, other):
        if isinstance(other, Vector):
            added = [ a + b for a, b in zip(self, other) ]
        elif isinstance(other, (int, float)):
            added = [ a + other for a in self ]
        else:
            raise ValueError("Addition with type {} not supported".format(type(other)))
        return self.__class__(added)

    def __truediv__(self, other):
        if isinstance(other, Vector):
            divided = [ self[i] / other[i] for i in range(len(self)) ]
        elif isinstance(other, (int, float)):
            divided = [ a / other for a in self ]
        else:
            raise ValueError("Division with type {} not supported".format(type(other)))

        return self.__class__(divided)

    def __sub__(self, other):
        """ Returns the vector difference of self and other """
        if isinstance(other, Vector):
            subbed = [ a - b for a, b in zip(self, other) ]
        elif isinstance(other, (int, float)):
            subbed = [ a - other for a in self ]
        else:
            raise ValueError("Subtraction with type {} not supported".format(type(other)))
        return self.__class__(subbed)

    def __mul__(self, other):
        if isinstance(other, Vector):
            product = sum([a * b for a, b in zip(self, other) ])
            return product
        elif isinstance(other, (int, float)):
            product = [ a * other for a in self ]
            return self.__class__(product)
        else:
            raise ValueError("Multiplication with type {} not supported".format(type(other)))

    def __iter__(self):
        return self.values.__iter__()
    def __len__(self):
        return len(self.values)
    def __getitem__(self, key):
        return self.values[key]
    def __repr__(self):
        return str(self.values)

def proj(a,v):
    if isinstance(a,Vector) and isinstance(v,Vector):
        return a*((a*v)/(a*a))
    else:
        raise ValueError("Projection with these types is not supported.")

def angle(a,b):
    if isinstance(a, Vector) and isinstance(b,Vector):
        res =  (a*b)/(a.norm()*b.norm())
        if res > 1 or res < -1:
            return arccos(sign(res))
        else:
            return arccos(res)
    else:
        raise ValueError("Angle with these types is not defined.")

def cross(a,b):
    if isinstance(a, Vector) and isinstance(b,Vector):
        return Vector([a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]])
    else:
        raise ValueError("Cross product with these types is not supported.")
