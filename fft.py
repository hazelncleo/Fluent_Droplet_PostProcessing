import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft, fftfreq, fftshift, fft2, rfft2, rfftfreq, rfft
from scipy.interpolate import LinearNDInterpolator
import skimage as ski
from photutils.profiles import RadialProfile
import glob as glob

def f_1(x):
    return 3*np.cos(2*np.pi*x)

def f_2(x):
    return 5*np.sin(50*np.pi*x/4)

def f_3(x):
    return 10*np.cos((2*np.pi/10)*x)

def F_1(X,Y):
    return np.cos(2*np.pi*X)

def F_2(X,Y):
    return np.sin(25*np.pi*Y/2)

def F_3(X,Y):
    return np.cos(np.pi*(X+2*Y))

def F_4(X,Y):
    return np.sin(2*np.pi*X/10)


def test_fft():
    samples_per = 25
    n = 20
    n_samples = (n*samples_per)+1
    sample_spacing = 1/samples_per
    x = np.linspace(0,n,n_samples)
    y = np.linspace(0,n,n_samples)
    X,Y = np.meshgrid(x,y)
    window = ski.filters.window(('tukey', 0.75), X.shape)

    z_1d = f_1(x) + f_2(x) + f_3(x) + np.random.normal(0,1,(n_samples))
    z_2d = F_1(X,Y) + F_2(X,Y) + F_3(X,Y) + F_4(X,Y) + np.random.normal(0,1,X.shape)

    zf = fft(z_1d)

    xf = fftfreq(n_samples, sample_spacing)
    yf = fftfreq(n_samples, sample_spacing)
    yf_1d = fftshift(yf)
    zf_1d = fftshift(zf)

    Xf, Yf = np.meshgrid(xf,yf)
    
    zf_2d_window = fftshift(fft2(z_2d*window))
    zf_2d_window = np.log(np.abs(zf_2d_window)**2)
    Xf, Yf = fftshift(Xf), fftshift(Yf)

    f,(ax_1,ax_2) = plt.subplots(2,4, width_ratios=[1,1,1,1], height_ratios=[1,1], figsize=(20,10))

    ax_1[0].plot(x, z_1d, '-k')

    ax_1[1].plot(yf_1d, np.abs(zf_1d)/np.max(np.abs(zf_1d)), '-r')

    ax_1[2].plot(1/yf[yf != 0],np.abs(zf[yf != 0])/np.max(np.abs(zf[yf != 0])), '-g')

    ax_2[0].pcolormesh(X,Y,z_2d*window,shading='gouraud')

    ax_2[1].pcolormesh(Xf,Yf,zf_2d_window,shading='gouraud')

    frequency= np.sqrt(np.power(Xf,2) + np.power(Yf,2)).flatten()
    truth_array = np.logical_and(frequency > 0.01 , frequency <= 0.6) 
    zf_2d_window = zf_2d_window.flatten()

    ax_2[2].scatter(1/frequency[truth_array], zf_2d_window[truth_array])

    plt.show()


def main():
    
    # Set sampling data
    elements_along_width = 32
    elements_along_length = 320
    sample_spacing = 500/elements_along_length
    N = elements_along_length + 1
    M = elements_along_width + 1

    files = glob.glob('.\\files\\*.npy')

    for i,file in enumerate(files):
        # Read isosurface coordinates
        with open(file, 'rb') as f:
            coords = np.load(f)*1e6

        # Centre Z coordinates around 0
        coords[:,2] = coords[:,2] - np.mean(coords[:,2])
        
        # Create plot
        f,ax = plt.subplots(2,4, width_ratios=[1,1,1,1], height_ratios=[1,1], figsize=(20,10))

        # Create Meshes and other helpers
        channel_width = np.linspace(-25,25,M)
        channel_length = np.linspace(-250,250,N)
        X_vert,Y_vert = np.meshgrid(channel_width,channel_length)
        X_hor,Y_hor = np.meshgrid(channel_length,channel_width)
        window = ski.filters.window(('tukey', 0.75), X_vert.shape)
        interp = LinearNDInterpolator(coords[:,:2], coords[:,2])

        # Interpolate data onto ordered grid
        Z_vert = interp(X_vert,Y_vert)
        Z_hor = interp(X_hor,Y_hor)

        # Get frequencies
        xf = fftfreq(M, sample_spacing)
        yf = fftfreq(N, sample_spacing)
        Xf,Yf = np.meshgrid(xf,yf)

        # 2D FFT with and without window functions
        Zf_vert_w = fft2(Z_vert*window)
        Zf_hor_w = fft2(Z_hor*np.transpose(window))

        Xf, Yf, Zf_vert_w, Zf_hor_w = fftshift(Xf), fftshift(Yf), fftshift(Zf_vert_w), fftshift(Zf_hor_w)

        # Get PSD
        Zf_vert_w = np.log(np.abs(Zf_vert_w)**2)
        Zf_hor_w = np.log(np.abs(Zf_hor_w)**2)

        truth_array = (coords[:,0] > -25.00001) & (coords[:,0] < 25.00001)
        ax[0,0].pcolormesh(X_vert,Y_vert,Z_vert,shading='gouraud')
        ax[0,0].set_xlim(-250,250)
        ax[0,0].set_ylim(-250,250)
        
        ax[0,1].pcolormesh(X_hor,Y_hor,Z_hor)
        ax[0,1].set_xlim(-250,250)
        ax[0,1].set_ylim(-250,250)

        ax[0,2].pcolormesh(X_vert,Y_vert,Z_vert*window)
        ax[0,2].set_xlim(-250,250)
        ax[0,2].set_ylim(-250,250)

        ax[0,3].pcolormesh(X_hor,Y_hor,Z_hor*np.transpose(window))
        ax[0,3].set_xlim(-250,250)
        ax[0,3].set_ylim(-250,250)

        ax[1,0].pcolormesh(Xf,Yf,Zf_vert_w)
        ax[1,0].set_xlim(-10*np.max(Xf),10*np.max(Xf))

        ax[1,1].pcolormesh(np.transpose(Yf),np.transpose(Xf),Zf_hor_w,shading='gouraud')
        ax[1,1].set_ylim(-10*np.max(Xf),10*np.max(Xf))

        frequency= np.sqrt(np.power(Xf,2) + np.power(Yf,2)).flatten()
        truth_array = np.logical_and(frequency > 0.022 , frequency <= 0.4) 
        Zf_vert_w = Zf_vert_w.flatten()
        ax[1,2].scatter(1/frequency[truth_array], Zf_vert_w[truth_array],s=0.75,c='k')
        for mult in [0.5,1,2,3]:
            ax[1,2].axvline(8.8*mult,0,1, color='r',lw=1,alpha=0.5,ls='--')
            ax[1,3].axvline(8.8*mult,0,1, color='r',lw=1,alpha=0.5,ls='--')

        frequency = np.transpose(np.sqrt(np.power(Xf,2) + np.power(Yf,2))).flatten()
        truth_array = np.logical_and(frequency > 0.022 , frequency <= 0.4) 
        Zf_hor_w = Zf_hor_w.flatten()
        ax[1,3].scatter(1/frequency[truth_array], Zf_hor_w[truth_array],s=0.75,c='k')

        plt.savefig('.\\images\\image_{}.png'.format(i),dpi=750)
    




main()