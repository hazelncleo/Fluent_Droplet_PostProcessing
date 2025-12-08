import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft, fftfreq, fftshift, fft2, rfft2, rfftfreq, rfft
from scipy.interpolate import LinearNDInterpolator


def f_1(x):
    return 3*np.cos(2*np.pi*x)

def f_2(x):
    return 5*np.sin(2*np.pi*x/4)

def f_3(x):
    return 10*np.cos(2*np.pi*x/10)

def F_1(X,Y):
    return np.cos(np.pi*X)

def F_2(X,Y):
    return np.sin(np.pi*Y)

def F_3(X,Y):
    return np.cos(np.pi*(X+Y))



def main():
    test = True
    if test:
        samples_per = 50
        n = 10
        n_samples = (n*samples_per)
        sample_spacing = 1/samples_per
        x = np.linspace(0,n,n_samples+1)[:-1]
        y = np.linspace(0,n,n_samples+1)[:-1]
        X,Y = np.meshgrid(x,y)

        z_1d = f_1(x) + f_2(x) + f_3(x) + np.random.normal(0,3,(n_samples))
        z_2d = F_3(X,Y) + F_1(X,Y) + F_2(X,Y) + np.random.normal(0,1,X.shape)# + 5*np.sin(2*np.pi*Y/4)# + 10*np.cos(np.pi*(X+Y)/6)

        zf_1d = rfft(z_1d)[1:]

        xf = fftfreq(n_samples, sample_spacing)
        yf = rfftfreq(n_samples, sample_spacing)
        yf_1d = yf[1:]

        f,(ax_1,ax_2) = plt.subplots(2,3)

        ax_1[0].plot(x, z_1d, '-k')
        ax_1[1].plot(yf_1d, np.abs(zf_1d)/np.max(np.abs(zf_1d)), '-r')
        ax_1[2].plot(1/yf_1d,np.abs(zf_1d)/np.max(np.abs(zf_1d)), '-g')
        
        yf = fftfreq(n_samples, sample_spacing)
        Xf,Yf = np.meshgrid(xf,yf)
        
        
        zf_2d = fft2(z_2d)
        zf_2d_sliced = fftshift(zf_2d)
        Xf_sliced = fftshift(Xf)
        Yf_sliced = fftshift(Yf)
        Xf_sliced = np.delete(Xf_sliced,int(samples_per*n/2),axis=0)
        Xf_sliced = np.delete(Xf_sliced,int(samples_per*n/2),axis=1)
        Yf_sliced = np.delete(Yf_sliced,int(samples_per*n/2),axis=0)
        Yf_sliced = np.delete(Yf_sliced,int(samples_per*n/2),axis=1)
        zf_2d_sliced = np.delete(zf_2d_sliced,int(samples_per*n/2),axis=0)
        zf_2d_sliced = np.delete(zf_2d_sliced,int(samples_per*n/2),axis=1)
        
        #normed_wl = np.sqrt(np.power(1/Xf,2) + np.power(1/Yf,2)).flatten()
        #flat_zf_2d = zf_2d.flatten()
        #norm_zf_2d = np.abs(flat_zf_2d)


        ax_2[0].imshow(z_2d)
        #ax_2[1].pcolormesh(fftshift(Xf),fftshift(Yf),np.abs(zf_2d)/np.max(np.abs(zf_2d)),vmin=0,vmax=1, shading='gouraud')
        #ax_2[1].imshow(np.abs(fftshift(zf_2d))/np.max(np.abs(zf_2d)),vmin=0,vmax=1)
        ax_2[1].pcolormesh(fftshift(Xf),fftshift(Yf),np.log(np.abs(fftshift(zf_2d))**2),shading='gouraud')
        ax_2[2].pcolormesh(1/Xf_sliced,1/Yf_sliced,np.log(np.abs(fftshift(zf_2d_sliced))**2),shading='gouraud')
        #ax_2[2].pcolormesh(1/Xf,1/Yf,np.abs(zf_2d)/np.max(np.abs(zf_2d)),vmin=0,vmax=1, shading='gouraud')
        #ax_2[2].scatter(normed_wl,norm_zf_2d/np.max(norm_zf_2d))
        
    

    else: 
        with open('outfile.npy', 'rb') as f:
            coords = np.load(f)
        N = 641
        M = 65
        sample_spacing = 50/64

        f,ax = plt.subplots(2,2)
        #ax[2].remove()
        #ax[2] = f.add_subplot(1, 3, 3, projection='3d')

        x = np.linspace(-25e-6,25e-6,M)
        y = np.linspace(-250e-6,250e-6,N)
        X,Y = np.meshgrid(x,y)
        
        interp = LinearNDInterpolator(coords[:,:2], coords[:,2])
        Z = interp(X,Y)

        xf = fftfreq(M, sample_spacing)*64/55
        yf = rfftfreq(N, sample_spacing)*64/55
        Xf,Yf = np.meshgrid(xf,yf)
        zf = rfft2(Z)[1:,1:]
        Xf = Xf[1:,1:]
        Yf = Yf[1:,1:]
        normed_wl = np.sqrt(np.power(1/Xf,2) + np.power(1/Yf,2)).flatten()
        flat_zf = zf.flatten()
        zf_norm = np.abs(flat_zf)

        zf_cropped = zf[9:,9:]
        Xf_cropped = Xf[9:,9:]
        Yf_cropped = Yf[9:,9:]
        normed_wl_cropped = np.sqrt(np.power(1/Xf_cropped,2) + np.power(1/Yf_cropped,2)).flatten()
        flat_zf_cropped = zf_cropped.flatten()
        zf_norm_cropped = np.abs(flat_zf_cropped)

        truth_array = (coords[:,0] > -25.00001e-6) & (coords[:,0] < 25.00001e-6)
        ax[0,0].scatter(coords[truth_array][:,0],coords[truth_array][:,1],s=0.75,c=coords[truth_array][:,2],vmin=np.min(Z),vmax=np.max(Z))
        ax[0,0].set_xlim(-250e-6,250e-6)
        ax[0,0].set_ylim(-250e-6,250e-6)
        ax[0,1].set_xlim(-250e-6,250e-6)
        ax[0,1].set_ylim(-250e-6,250e-6)

        ax[0,1].scatter(X,Y,s=0.75,c=Z,vmin=np.min(Z),vmax=np.max(Z))
        #levels = np.linspace(0, 1, 7)
        #ax[2].scatter(Xf,Yf,s=0.75,c=(1/np.max(zf_norm))*zf_norm)
        ax[1,0].scatter(normed_wl,zf_norm/np.max(zf_norm))
        ax[1,1].scatter(normed_wl_cropped, zf_norm_cropped/np.max(zf_norm_cropped))
        #ax[2].plot_surface(Xf,Yf,zf_norm)#,levels=levels)
    





    plt.show()
    




main()