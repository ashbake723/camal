# Plotting Functions

plt.plot(time[1:][::3],flux[1:][::3]/texp[1:][::3],'ro')
plt.plot(time[2:][::3],flux[2:][::3]/texp[2:][::3],'go')
plt.plot(time[::3],flux[::3]/texp[::3],'bo')



for i in range(len(flux[::3])):
     if all(flux[j:j+3]>0):
         ts, tdf = getpt(flux[j],flux[j+1],flux[j+2])
         slope.append(ts)
         df.append(tdf)
     else:
         slope.append(0)
         df.append(0)
     j+=3



out = signal.fftconvolve(camal.gaussian(np.arange(len(df)),len(df)/2,10),df,'same')
