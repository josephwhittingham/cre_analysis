import numpy as np
import sys

class logspace:
	""" Class for functions which are best represented as piecewise powerlaws """	

	def __init__(self, x, y, axis=None, extrapolate=False, ctype=np.float64, errstate='ignore'):
		"""
	    Calculate piecewise power laws necessary for interpolation (+ extrapolation) and integration.
		The x array is interpreted as bin boundaries between which a piece-wise power law functions are calculated.
		

		Parameters
		----------
		x : numpy array
			The 'x' points at which 'y' is sampled.

		y : numpy array
			Array to be integrated

		axis : int, optional
			Axis along which the piece-wise power laws are calculated. Default is last axis.

		ctype : type, optional
			Type to be used for calculation. Default is np.float64.

		errstate : {'ignore', 'warn', 'raise', 'call', 'print', 'log'}, optional
			Set how floating-point errors are handled. Default is 'ignore'.

		Return
		------
		numpy array
			Result of integration
		"""

		self.ctype = ctype

		# Select the correct axis
		if axis == None:
			self.axis = len(y.shape) -1
		elif type(axis) is not int:
			sys.exit("Parameter choice 'axis = {:}' is not allowed.".format(axis))
		else:
			self.axis = axis

		# Convert x array to numpy type
		if type(x) is int or type is float:
			self.x = np.array([x])
		elif type(x) is list:
			self.x = np.array(x)
		else:
			self.x = np.copy(x)

		if len(self.x.shape) > 1:
			print("Warning: Interpolation, extrapolation, and integration to arbitrary boundaries works only for 1D x arrays")

		self._x_orig = np.copy(x) # Backup copy necessary for interpolation/extrapolation

		# Convert y array to numpy type
		if type(y) is int or type is float:
			self.y = np.array([y])
		elif type(y) is list:
			self.y = np.array(y)
		else:
			self.y = np.copy(y)

		# Expand dimensions of x array
		for i in np.arange(len(self.x.shape), len(self.y.shape) - self.axis):
			self.x = np.expand_dims(self.x, i)

		self.x = np.broadcast_to(self.x, self.y.shape)

		# Check if x array is sorted
		if np.any(np.diff(self.x, axis=self.axis) < 0):
			sys.exit("Array x is not sorted!")

		if self.y.shape[self.axis] <= 1:
			return 0

		with np.errstate(all=errstate):
			self.alpha = np.log(self.__lval(self.y) / self.__rval(self.y)) / np.log( self.__lval(self.x) / self.__rval(self.x)).astype(self.ctype)
			self.C = self.__lval(self.y) * np.power(self.__lval(self.x), -self.alpha)
			self.nonfinite = np.where(np.logical_not(np.logical_and(np.isfinite(self.alpha), np.isfinite(self.C))))
			self.alpha[self.nonfinite] = 0
			self.C[self.nonfinite] = 0

	def __lval(self, arr):
		""" Returns values of left bin boundaries """
		return np.take(arr, np.arange(arr.shape[self.axis] - 1), axis=self.axis)

	def __rval(self, arr):
		""" Returns values of left bin boundaries """
		return np.take(arr, np.arange(1, arr.shape[self.axis]), axis=self.axis)


	def integrate(self, x_min = None, x_max = None, extrapolate=False, rtype=np.float64):
		"""
		Integrate y(x) using samples with piece-wise power
		laws. Integration ranges around zero values in y are treated as
		zero.

		Parameters
		----------
		x_min : float, optional
		    Lower integration boundary. Set extrapolate=True if value
		    outside of original 'x' array should be considered.

		x_max : float, optional
		    Upper integration boundary. Set extrapolate=True if value
		    outside of original 'x' array should be considered.

		rtype : type, optional
			Type to be used for return. Default is np.float64.

		Return
		------
		numpy array
			Result of integration
		"""

		if x_min is None and x_max is None:
			return np.sum(np.where(self.alpha!= -1,
								   self.C / (self.alpha + 1) * ( np.power(self.__rval(self.x), self.alpha + 1) - np.power(self.__lval(self.x), self.alpha+1)),
								   self.C * np.log(self.__rval(self.x)/self.__lval(self.x))), axis=self.axis).astype(rtype)
		else:
			# Works only for 1D x_int arrays
			if x_min == None:
				x_min = np.array([self._x_orig[0]])
			elif type(x_min) is list or type(x_min) is np.ndarray:
				if len(x_min) < 1 or len(x_min) > 1:
					sys.exit("Only scalar or array of 1D array of length 1 is allowed for x_min")
			else:
				x_min = np.array([x_min])

			if x_max == None:
				x_max = np.array([self._x_orig[-1]])
			elif type(x_max) is list or type(x_max) is np.ndarray:
				if len(x_max) < 1 or len(x_max) > 1:
					sys.exit("Only scalar or array of 1D array of length 1 is allowed for x_max")
			else:
				x_max = np.array([x_max])

			if x_max < x_min:
				sys.exit("x_max is smaller than x_min")

			# Define an interpolating x array including new boundaries
			x = np.concatenate((x_min, self._x_orig[np.logical_and( self._x_orig > x_min, self._x_orig < x_max)], x_max))

			# Only works for 1D x and original x array
			bin_index = np.digitize(x, self._x_orig) # bin index for every interpolation x point

			if extrapolate:
				# Values lying outside of original x range are now mapped to outer bins
				bin_index[np.where(bin_index == 0)] = 1
				bin_index[np.where(bin_index == len(self._x_orig))] = len(self._x_orig) - 1
			else:
				# Values lying on rightmost bin boundary are included
				bin_index[np.where(x == self._x_orig[-1])] = len(self._x_orig) - 1

			# Indices of x array and of bins inside the integration range
			inside_range = np.where(np.logical_and(bin_index > 0, bin_index < len(self._x_orig)))[0]
			bins_inside = bin_index[inside_range] - 1

			x_shape = list(self.y.shape)
			x_shape[self.axis] = len(x)

			for i in np.arange(len(x.shape), len(self.y.shape) - self.axis):
				x = np.expand_dims(x, i)

			x = np.broadcast_to(x, x_shape)

			x = np.take(x, inside_range, axis=self.axis)
			C = np.take(self.C, bins_inside[:-1], axis=self.axis)
			alpha = np.take(self.alpha, bins_inside[:-1], axis=self.axis)

			return np.sum(np.where(alpha!= -1,
								   C / (alpha + 1) * ( np.power(self.__rval(x), alpha + 1) - np.power(self.__lval(x), alpha+1)),
								   C * np.log(self.__rval(x)/self.__lval(x))), axis=self.axis).astype(rtype)


	def interpolate(self, x_int, extrapolate=False, rtype=np.float64):
		"""
		Calculate an interpolation array

		Parameters
		----------
		x_int : numpy array
			The points at which interpolated values are calculated.

		rtype : type, optional
			Type to be used for return. Default is np.float64.

		Return
		------
		numpy array
			Result of integration
		"""
		
		scalar_input = False
		
		if type(x_int) is int or type(x_int) is float:
			x_int = np.array([x_int])
			scalar_input = True
		elif type(x_int) is list:
			x_int = np.array(x_int)

		if len(x_int.shape) > 1:
			sys.exit("Only 1D arrays are allowed for 'x_int'")

		# Only works for 1D x_int and original x array
		bin_index = np.digitize(x_int, self._x_orig) # bin index for every interpolation x point


		if extrapolate:
			bin_index[np.where(bin_index == 0)] = 1
			bin_index[np.where(bin_index == len(self._x_orig))] = len(self._x_orig) - 1
		else:
			# Values lying on rightmost bin boundary are included
			bin_index[np.where(x_int == self._x_orig[-1])] = len(self._x_orig) - 1

		inside_range = np.where(np.logical_and(bin_index > 0, bin_index < len(self._x_orig)))[0]
		bins_inside = bin_index[inside_range] - 1

		res_shape = list(self.y.shape)
		res_shape[self.axis] = len(x_int)
		res = np.zeros(res_shape, dtype=self.ctype)

		if len(self.y.shape) == 1:
			res[inside_range] = self.C[bins_inside] * np.power(x_int[inside_range], self.alpha[bins_inside])

			if scalar_input:
				# Return scalar if interpolation point is given as scalar
				res = res[0]

		elif len(self.y.shape) == 2:
			if self.axis == 0:
				res[inside_range, :] = self.C[bins_inside, :] * np.power(x_int[inside_range, np.newaxis], self.alpha[bins_inside, :])
			else:
				res[:, inside_range] = self.C[:, bins_inside] * np.power(x_int[inside_range], self.alpha[:, bins_inside])

		elif len(self.y.shape) == 3:
			if self.axis == 0:
				res[inside_range, :, :] = self.C[bins_inside, :, :] * np.power(x_int[inside_range, np.newaxis,  np.newaxis], self.alpha[bins_inside, :, :])
			elif self.axis == 1:
				res[:, inside_range, :] = self.C[:, bins_inside, :] * np.power(x_int[np.newaxis, inside_range, np.newaxis], self.alpha[:, bins_inside, :])
			else:
				res[:, :, inside_range] = self.C[:, :, bins_inside] * np.power(x_int[np.newaxis, np.newaxis, inside_range], self.alpha[:, :, bins_inside])
		else:
			sys.exit("Only up to 3 dimensions are supported for y")
		
		return res
	
		



def logtrapz(y, x, axis=None, ctype=np.float128, rtype=np.float64, errstate='ignore'):
	"""
	Integrate y(x) using samples with piece-wise power
	laws. Integration ranges around zero values in y are treated as
	zero.

	Parameters
	----------
	y : numpy array
	    Array to be integrated

	x : numpy array
	    The 'x' points at which 'y' is sampled.

	axis : int, optional
	    Axis along which to integrate. Default is the last axis

	ctype : type, optional
	    Type to be used for calculation. Default is np.float128.

	rtype : type, optional
	    Type to be used for return. Default is np.float64.

	errstate : {'ignore', 'warn', 'raise', 'call', 'print', 'log'}, optional
	    Set how floating-point errors are handled. Default is 'ignore'.
	    

	Return
	------
	numpy array
	    Result of integration
	"""
	
	if axis == None:
		axis = len(y.shape) -1
	elif type(axis) is not int:
		sys.exit("Parameter choice 'axis = {:}' is not allowed.".format(axis))

	if type(x) is int or type is float:
		x = np.array([x])
	elif type(x) is list:
		x = np.array(x)

	for i in np.arange(len(x.shape), len(y.shape) - axis):
		x = np.expand_dims(x, i)

	x = np.broadcast_to(x, y.shape)

	if y.shape[axis] <= 1:
		return 0

	def lval(arr, axis=axis):
		return np.take(arr, np.arange(arr.shape[axis] - 1), axis=axis)

	def rval(arr, axis=axis):
		return np.take(arr, np.arange(1, arr.shape[axis]), axis=axis)

	with np.errstate(all=errstate):
		alpha = np.log(lval(y, axis=axis) / rval(y, axis=axis)) / np.log( lval(x, axis=axis) / rval(x, axis=axis)).astype(ctype)
		C = lval(y, axis=axis) * np.power(lval(x, axis=axis), -alpha)
		nonfinite = np.where(~(np.isfinite(alpha) & np.isfinite(C)))
		alpha[nonfinite] = 0
		C[nonfinite] = 0
			
	return np.sum(np.where(alpha!= -1,
						   C / (alpha + 1) * ( np.power(rval(x, axis=axis), alpha + 1) - np.power(lval(x, axis=axis), alpha+1)),
						   C * np.log(rval(x)/lval(x))), axis=axis).astype(rtype)


def loginterp(y, x, x_int, axis=None, extrapolate=False, ctype=np.float128, rtype=np.float64, errstate='ignore'):
	"""
	Calculate an interpolation array

	Parameters
	----------
	y : numpy array
	    Array to be interpolated.

	x : numpy array
	    The 'x' points at which 'y' is sampled.

	x_int : numpy array
	    The points at which interpolated values are calculated.

	axis : int, optional
	    Axis along which to interpolate. Default is the last axis

	ctype : type, optional
	    Type to be used for calculation. Default is np.float128.

	rtype : type, optional
	    Type to be used for return. Default is np.float64.

	errstate : {'ignore', 'warn', 'raise', 'call', 'print', 'log'}, optional
	    Set how floating-point errors are handled. Default is 'ignore'.
	    

	Return
	------
	numpy array
	    Result of integration
	"""

	if type(x) is int or type(x) is float:
		x = np.array([x])
	elif type(x) is list:
		x = np.array(x)

	if type(x_int) is int or type(x_int) is float:
		x_int = np.array([x_int])
	elif type(x_int) is list:
		x_int = np.array(x_int)

	if len(x_int.shape) > 1:
		sys.exit("Only 1D arrays are allowed for 'x_int'")

	if len(x.shape) > 1:
		sys.exit("Only 1D arrays are allowed for 'x'")
	
	if axis == None:
		axis = len(y.shape) -1
	elif type(axis) is not int:
		sys.exit("Parameter choice 'axis = {:}' is not allowed.".format(axis))

	x_orig = np.copy(x)

	for i in np.arange(len(x.shape), len(y.shape) - axis):
		x = np.expand_dims(x, i)

	x = np.broadcast_to(x, y.shape)		
	def lval(arr, axis=axis):
		return np.take(arr, np.arange(arr.shape[axis] - 1), axis=axis)

	def rval(arr, axis=axis):
		return np.take(arr, np.arange(1, arr.shape[axis]), axis=axis)

	with np.errstate(all=errstate):
		alpha = np.log(lval(y, axis=axis) / rval(y, axis=axis)) / np.log( lval(x, axis=axis) / rval(x, axis=axis)).astype(ctype)
		C = lval(y, axis=axis) * np.power(lval(x, axis=axis), -alpha)
		nonfinite = np.where(~(np.isfinite(alpha) & np.isfinite(C)))
		alpha[nonfinite] = 0
		C[nonfinite] = 0

	# Only works for 1D x_int and original x array
	bin_index = np.digitize(x_int, x_orig) # bin index for every gas cell
	
	if extrapolate:
		bin_index[np.where(bin_index == 0)] = 1
		bin_index[np.where(bin_index == len(x_orig))] = len(x_orig) - 1
	else:
		# Values lying on rightmost bin boundary are included
		bin_index[np.where(x_int == x_orig[-1])] = len(x_orig) - 1

		
	inside_range = np.where(np.logical_and(bin_index > 0, bin_index < len(x_orig)))[0]
	bins_inside = bin_index[inside_range] - 1
		
	res_shape = list(y.shape)
	res_shape[axis] = len(x_int)
	res = np.zeros(res_shape, dtype=ctype)

	if len(y.shape) == 1:
		res[inside_range] = C[bins_inside] * np.power(x_int[inside_range], alpha[bins_inside])

	elif len(y.shape) == 2:
		if axis == 0:
			res[inside_range, :] = C[bins_inside, :] * np.power(x_int[inside_range, np.newaxis], alpha[bins_inside, :])
		else:
			res[:, inside_range] = C[:, bins_inside] * np.power(x_int[inside_range], alpha[:, bins_inside])

	elif len(y.shape) == 3:
		if axis == 0:
			res[inside_range, :, :] = C[bins_inside, :, :] * np.power(x_int[inside_range, np.newaxis,  np.newaxis], alpha[bins_inside, :, :])
		elif axis == 1:
			res[:, inside_range, :] = C[:, bins_inside, :] * np.power(x_int[np.newaxis, inside_range, np.newaxis], alpha[:, bins_inside, :])
		else:
			res[:, :, inside_range] = C[:, :, bins_inside] * np.power(x_int[np.newaxis, np.newaxis, inside_range], alpha[:, :, bins_inside])
	else:
		sys.exit("Only up to 3 dimensions are supported for y")


	return res
