# Get the centroid position of the star

class get_centroid:
    def __init__(self):
        self.info = 'sourceFinder(image)'
 
	def sourceFinder(self,image):
	    ' Input: loaded data, Output: x,y center of aperture '
    	# Convolve image with gaussian kernel
    	kernel = np.outer(signal.gaussian(50,8), signal.gaussian(50,8))
    	blurred = signal.fftconvolve(image, kernel, mode='same')
   
    	# Take the normalized STD along x,y axes
    	xstd = np.std(blurred,axis=0)
    	ystd = np.std(blurred,axis=1)
    	xstdn = (xstd - np.median(xstd[200:300]))/max(xstd)
    	ystdn = (ystd - np.median(ystd[200:300]))/max(ystd)
    
    	# Determine center by maximum. Eventually add check that there's only one source!
    	try: x,y = np.where(xstdn == max(xstdn))[0][0], np.where(ystdn == max(ystdn))[0][0]
    	except IndexError:
        	x,y = 0,0

    	self.x = x
        self.y = y