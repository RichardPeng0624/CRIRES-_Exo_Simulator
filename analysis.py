#The class focuses on processing theoretical spectra of planets.
import glob
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import fftpack,signal
import matplotlib
from matplotlib.ticker import ScalarFormatter
from astropy.visualization import (MinMaxInterval, SqrtStretch, ImageNormalize, ZScaleInterval)
import scipy
import json
import subprocess
import os
  
#The class focuses on processing theoretical spectra of planets.  

class planet_para:
    '''
    #FORMAT SUPPORT:     
    #disper_planet, disper_star, disper_sky: pd.DataFrame
    #data_planet, data_star, data_sky: np.ndarray or pd.DataFrame
    '''
    def __init__(self, disper_planet=None, disper_star=None,\
                 data_planet=None, data_star=None, disper_sky=None, data_sky=None):
        
        # disper_planet, disper_star, disper_sky: the dataset of the 2-D flux distribution of the planet, the star and the sky, in pd.Dataframe
        # data_planet, data_star, data_sky: the dataset of the 1-D extracted spectrum of the star, the planet, and the sky
        
        self.disper_planet = disper_planet
        self.data_planet = data_planet
        
        self.disper_star = disper_star
        self.data_star = data_star
        
        self.disper_sky = disper_sky
        self.data_sky = data_sky
        

        
    #-------------------
    #This function aims to calculate the color indexes of two given band for the spectrum.
    
    def color_index(self,\
        band1_w=None, band2_w=None, transmission1=None, transmission2=None,\
        zero1=None, zero2=None, atmosphere='False',\
        c1=None, c2=None,\
        airmass=None):
        
        '''
        #This function aims to calculate the color indexes of two given band for the spectrum.
        # band1_w, band2_w: two selected bands 
        # transmission1, transmission2: transmission profiles for the two bands
        # zero1, zero2: zero flux for two bands. ATTENZIONE: flux not magnitude!
        
        # --Atmospheric extinction parameters
        # --- atmosphere: The default setting is 'False', means no atmospheric extinction is considered. 
        # --- c1, c2: atmospheric extinction coefficients
        # --- airmass: only need to be defined when you would like to calculate the atmospheric extinction.
        '''
        
        if isinstance(self.data_planet, np.ndarray):

            wave = self.data_planet[0]
            flux = self.data_planet[1]

            mask_1=np.array((wave>band1_w[0])&(wave<band1_w[1]))
            mask_2=np.array((wave>band2_w[0])&(wave<band2_w[1]))

            data_1=np.array([self.data_planet[0][mask_1],self.data_planet[1][mask_1]])
            transmission_1=np.interp(data_1[0], transmission1[0], transmission1[1])
            
            data_2=np.array([self.data_planet[0][mask_2],self.data_planet[1][mask_2]])
            transmission_2=np.interp(data_2[0], transmission2[0], transmission2[1])

            int_flux_1= np.trapz(data_1[1]*transmission_1, data_1[0])/np.trapz(transmission_1, data_1[0])
            int_flux_2= np.trapz(data_2[1]*transmission_2, data_2[0])/np.trapz(transmission_2, data_2[0])
            
            

        elif isinstance(self.data_planet, pd.DataFrame):

            wave = self.data_planet['wave']
            flux = self.data_planet['flux']

            mask_1=np.array((wave>band1_w[0])&(wave<band1_w[1]))
            mask_2=np.array((wave>band2_w[0])&(wave<band2_w[1]))

            data_1=self.data_planet[mask_1]
            transmission_1=np.interp(data_1[0], transmission1[0], transmission1[1])
            
            data_2=self.data_planet[mask_2]
            transmission_2=np.interp(data_2[0], transmission2[0], transmission2[1])

            int_flux_1= np.trapz(data_1.flux*transmission_1, data_1.wave)/np.trapz(transmission_1,  data_1.wave)
            int_flux_2= np.trapz(data_2.flux*transmission_2, data_2.wave)/np.trapz(transmission_2,  data_2.wave)
        
        else:
            
            print ('Sweet but dumb, we do not support this format. Choose one from np.ndarray or pd.DataFrame.')
            
            
        #lam_eff_1=np.trapz(data_1.flux*transmission_1*(data_1.wave**2), data_1.wave)/np.trapz(data_1.flux*transmission_1*(data_1.wave), data_1.wave)
        #lam_eff_2=np.trapz(data_2.flux*transmission_2*(data_2.wave**2), data_2.wave)/np.trapz(data_2.flux*transmission_2*(data_2.wave), data_2.wave)
        
        c=scipy.constants.speed_of_light
        
        #zero_flux_1=zero1*1e-23*(c/(lam_eff_1**2))
        #zero_flux_2=zero2*1e-23*(c/(lam_eff_2**2))
        
        if atmosphere=='True':
            mag_1=-2.5*np.log10(int_flux_1/zero1)-c1*(airmass-1)
            mag_2=-2.5*np.log10(int_flux_2/zero2)-c2*(airmass-1)
        
        else:
        
            mag_1=-2.5*np.log10(int_flux_1/zero1)
            mag_2=-2.5*np.log10(int_flux_2/zero2)
              
        
        c_index=mag_1-mag_2
        #c_index=-2.5*np.log10(int_flux_1/int_flux_2)+(zero1-zero2)+(c1-c2)*(airmass-1.0)
        
        return (c_index, mag_1, mag_2, data_1, transmission_1, data_2, transmission_2)
    
    
    #----------------------
    def p_2_s (self,  order_num, plot=True, planet_posi= None, sky=False):
        '''
        This function aims to calculate both the intrinsic and the in-situ planet-to-star flux ratio 
        #if planet position is setted as False, only intrinsic p2s ratio is calculated 
        '''
        font = {'size': 4}
        plt.rcParams.update({'font.size': font['size']})
        
        if sky==True:
            self.disper_star+=self.disper_sky
            

        p_2_s_intr=np.float64(self.disper_planet['dat_diff_0']/self.disper_star['dat_diff_0'])
        
        p_2_s_intr=np.log10(p_2_s_intr)
        p_2_s_intr[np.isinf(p_2_s_intr)]=np.nan
        
        
        
        if planet_posi != False:
            
            p_2_s_planet = self.disper_planet['dat_diff_0']/self.disper_star['dat_diff_%s'%planet_posi]
            p_2_s_planet = np.log10(p_2_s_planet)
            p_2_s_planet[np.isinf(p_2_s_planet)]=np.nan
            
        else:
            p_2_s_planet = []
        
        
        
        if plot==True:

            fig,axs=plt.subplots(len(order_num),3)

            for i in range(len(order_num)):
                for j in range(3):
                    ax=axs[i][j]
                    ind=3*i+j
                    wave_range = self.disper_planet.loc[2045*ind:2045*(ind+1)]['wavelength(nm)']
                    mask_in=np.where(p_2_s_intr[2046*ind:2046*(ind+1)]==np.nan, True, False)
            
                    ax.plot(p_2_s_intr[2046*ind:2046*(ind+1)],\
                            linewidth=0.6, alpha=0.7)
                

                    if ind + 3 >= len(order_num)*3:
                        ax.set_xlabel('wavelength (nm)')
                    if j%3 ==0: 
                        ax.set_ylabel('$log_{10}(f_p/f_s)$')
                    if ind%3 !=0 :
                        ax.set_yticks([])
                        
                    ax.set_title('order:%s, detector:%s'%(order_num[i], j+1))
                    #ax.set_yscale('log')
                    
                    mi=np.nanmin(p_2_s_intr[2046*ind:2046*(ind+1)])
                    mi=np.round(mi,2)
                    ma=np.nanmax(p_2_s_intr[2046*ind:2046*(ind+1)])
                    ma=np.round(ma,2)
                    
                    ax.set_ylim(mi, ma*1.05)
                    ax.set_yticks(np.linspace(mi, ma, 3))
                    
                    tick=np.linspace(wave_range.min(),wave_range.max(), 3, dtype='int')
                    ax.set_xticks(ticks=np.linspace(2046*ind, 2046*(ind+1),3), labels=tick)
                    #ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))

            fig.suptitle('planet-to-star flux ratio')
            fig.subplots_adjust(hspace=1.3,wspace=0.3)
            plt.tight_layout()
            plt.close()
            
            
            if planet_posi != False:
                
                fig1,axs1=plt.subplots(len(order_num),3)

                for i in range(len(order_num)):
                    for j in range(3):
                        ax1=axs1[i][j]
                        ind=3*i+j
                        wave_range = self.disper_planet.loc[2045*ind:2045*(ind+1)]['wavelength(nm)']
                        mask_pl=np.where(p_2_s_planet[2046*ind:2046*(ind+1)]==np.nan, True, False)
                        
                        ax1.plot(p_2_s_planet[2046*ind:2046*(ind+1)],\
                                 linewidth=0.6, alpha=0.7)                     

                        if ind + 3 >= len(order_num)*3:
                            ax1.set_xlabel('wavelength (nm)')
                        if j%3 ==0: 
                            ax1.set_ylabel('$log_{10}(f_p/f_s)$') 
                            
                            
                        ax1.set_title('order:%s, detector:%s'%(order_num[i], j+1))
                        #ax1.set_yscale('log')
                                  
                        mi=np.nanmin(p_2_s_planet[2046*ind:2046*(ind+1)])
                        mi=np.round(mi,2)
                        ma=np.nanmax(p_2_s_planet[2046*ind:2046*(ind+1)])
                        ma=np.round(ma,2)
                        
                        ax1.set_ylim(mi, ma*1.05)
                        ax1.set_yticks(np.linspace(mi, ma, 3))
                        
                        tick=np.linspace(wave_range.min(),wave_range.max(),3, dtype='int')
                        ax1.set_xticks(ticks=np.linspace(2046*ind, 2046*(ind+1),3),labels=tick)
                        #ax1.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
                        
                        
                fig1.suptitle('planet-to-star flux ratio at planet position')
                fig1.subplots_adjust(hspace=1.3, wspace=0.3)
                plt.tight_layout()
                plt.close()
               
                
            
                return(p_2_s_intr, p_2_s_planet, fig, fig1)
            
            else:
        
                return (p_2_s_intr, p_2_s_planet, fig)
        
        else:
            
            return (p_2_s_intr, p_2_s_planet)
        
    #-------------------
    
    def mag_convert(self, flux_ratio, stellar_mag, error_b=0., error_p=0.):
        
        '''
        This function aims to convert the apparent magnitude of the planet based on the flux contrast ratio measured from 
        previous observations
        # error_b, error_p: the bottom and upper limit of error bar
        '''
        
        error=np.array([error_b, error_p])
        planet_mag_p = stellar_mag - 2.5*np.log10(flux_ratio+error[1])
        planet_mag_b = stellar_mag - 2.5*np.log10(flux_ratio-error[0])
        planet_mag = stellar_mag - 2.5*np.log10(flux_ratio)

        return ([planet_mag_b, planet_mag, planet_mag_p])
    
    #-------------------
    
    def mag_direct(self, radius, distance, band=None, transmission=None, zero=None, atmosphere='False', c=None, airmass=None):
                   
        '''
        --Convert the input paramters into CGS unit--

        # radius, distance: parameters of the planet

        # zero_point_flux: depending on the band and the photometric system, the suggestion is MKO potometric system.

        # band: the band for computation of magnitude.

        # spec_path: the path of input (theoretical) spectrum of the planet
        '''

        if isinstance(self.data_planet, np.ndarray):

            wave = self.data_planet[0]
            flux = self.data_planet[1]

            mask=np.array((wave>band[0])&(wave<band[1]))

            data_cut=np.array([wave[mask],flux[mask]])
            transmission=np.interp(data_cut[0], transmission[0], transmission[1])

            int_flux= np.trapz(data_cut[1]*transmission, data_cut[0])/np.trapz(transmission, data_cut[0])



        elif isinstance(self.data_planet, pd.DataFrame):

            wave = self.data_planet['wave']
            flux = self.data_planet['flux']

            mask=np.array((wave>band[0])&(wave<band[1]))

            data_cut=self.data_planet[mask]
            transmission=np.interp(data_cut[0], transmission[0], transmission[1])

            int_flux= np.trapz(data_cut.flux*transmission, data_cut.wave)/np.trapz(transmission,  data_cut.wave)

        else:

            print ('Sweet but dumb, we do not support this format. Choose one from np.ndarray or pd.DataFrame.')


        c=scipy.constants.speed_of_light
        pi=scipy.constants.pi

        dist_ratio=(radius/distance)**2

        F_rec=int_flux*dist_ratio


        if atmosphere=='True':
            mag=-2.5*np.log10(F_rec/zero)-c*(airmass-1)

        else:

            mag=-2.5*np.log10(F_rec/zero)



        return (mag, int_flux, F_rec)    
                   
     
