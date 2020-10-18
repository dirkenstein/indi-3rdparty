#include "nsdownload.h"
#include "kaf_constants.h"
#include "kai10100_constants.h"
#include <string.h>
#include <errno.h>
#include <stdlib.h>
#include  <unistd.h>
#include <string.h>
#include "nsdebug.h"
#include <math.h>
#include "indiccd.h"

void NsDownload::setFrameYBinning(int binning) {
			ctx->imgp->ybinning = binning;	

}

void NsDownload::setFrameXBinning(int binning) {
			ctx->imgp->xbinning = binning;	

}

void NsDownload::setImgSize(int siz) {
	rd->imgsz = siz;
}

void NsDownload::setNumExp(int n){
	ctx->nexp = n;		 	
}

void NsDownload::nextImage() {
		ctx->imgseq++;	
}
int NsDownload::getImgSeq(void) {
	return 		 ctx->imgseq;
	
}
void NsDownload::setSetTemp (float temp) {
		ctx->imgp->settemp = temp;	

}

void NsDownload::setActTemp(float temp) {
			ctx->imgp->acttemp = temp;	

}

void NsDownload::setExpDur(float exp){
				ctx->imgp->exp = exp;	
}

void NsDownload::setFrameType(int ft) {
	ctx->imgp->frame_type = ft;
}

void NsDownload::setIncrement(int inc) {
	ctx->increment = inc;	
}

void NsDownload::setFbase(const char * name) {
	strncpy(ctx->fbase,name, 64) ;	
}

void NsDownload::doDownload() {
	
	  ctx->imgp->expdate = time(NULL);
	  std::unique_lock<std::mutex> ulock(mutx);

		do_download = 1;
		go_download.notify_all();
}

bool NsDownload::inDownload() {
	return in_download||do_download;	
} 


bool NsDownload::getDoDownload() {
	return do_download;	
} 

unsigned char * NsDownload::getBuf() {
	if (retrBuf == NULL) return NULL;
	return retrBuf->buffer;	
}

int NsDownload::getBufSize() {
	if (retrBuf == NULL || retrBuf->buffer == NULL) {
		return 0;
	}
	return retrBuf->bufsiz;
}

size_t NsDownload::getBufImageSize() {
	if (retrBuf == NULL) return 0;

	return retrBuf->imgsz;	
}


void NsDownload::setImgWrite(bool w) {
	write_it = w;	
}

void NsDownload::freeBuf() {
	if (!retrBuf) return;
	if (retrBuf->buffer) free(retrBuf->buffer);
	retrBuf->buffer = NULL;
	retrBuf = NULL;
}

void NsDownload::setInterrupted(){
	std::unique_lock<std::mutex> ulock(mutx);
	interrupted = 1;
	go_download.notify_all();
}

void NsDownload::setZeroReads(int zeroes){
	zero_reads = zeroes;
}
 int NsDownload::getActWriteLines(){
		 return writelines;	
}

int NsDownload::get_max_x() {
	if (cn->getCamType() == kaf8300) {
		return KAF8300_MAX_X;
	} else if (cn->getCamType() == kai10100) {
		int res = KAI10100_MAX_X / ctx->imgp->xbinning;
		if (ctx->imgp->xbinning == 4) {
			res -= 1;
		}
		return res;
	} else {
		DO_ERR("unknown ccd type: %d", cn->getCamType());
		return 0;
	}
}

int NsDownload::downloader() 
{
			//struct ftdi_context * ftdid = cn->->getDataChannel();
			//struct ftdi_transfer_control * ctl;

      int rc2;
      int hardloop = 20; 
      int sleepage = 1000;
    
			int download =1;
			if (rd->nread > rd->bufsiz) {
            DO_ERR("image too large %d\n", rd->nread);
		     		return (-1);
			}
    	while((rc2 = cn->readData(rd->buffer+rd->nread, cn->getMaxXfer())) == 0 && hardloop > 0) {
				if (hardloop % 5  == 0) DO_INFO("W%d\n",hardloop);
    		usleep(sleepage);
				sleepage *= 2;
				if (sleepage > 100000) sleepage = 100000;
				hardloop--;
			}
      /* ctl = ftdi_read_data_submit(ftdid, rd->buffer+rd->nread,  maxxfer);
      if (ctl == NULL) {
      	DO_ERR( "unable to submit read: %d (%s)\n", rc2, ftdi_get_error_string(ftdid));
				return (-1);
      }
      DO_INFO( ".");
      rc2 = ftdi_transfer_data_done(ctl);
      */
   		if (rc2 < 0 ) {
        DO_ERR("unable to read download data: %d\n", rc2);
				return (-1);
			}
			rd->nread += rc2;
			if (rc2 != cn->getMaxXfer()) {
				DO_INFO("short! %d %d\n", rd->nblks, rc2);
			}		
			
			if (rc2 == 0) {
				readdone = 1;
			} else { rd->nblks ++; }
			if (rd->nread >= rd->imgsz) {
				readdone = 1;	
			}
			if (readdone) {
			  download=0;
				lastread = rc2;
			
				rb = rdd;
				retrBuf = &rb;
				rd->buffer = NULL;
			}	
			return download;		
}


/*

SIMPLE  =                    T                                                  
BITPIX  =                   16 /8 unsigned int, 16 & 32 int, -32 & -64 real     
NAXIS   =                    2 /number of axes                                  NAXIS1  =                 3326 /fastest changing axis                           
NAXIS2  =                 2504 /next to fastest changing axis                   BSCALE  =   1.0000000000000000 /physical = BZERO + BSCALE*array_value           
BZERO   =   32768.000000000000 /physical = BZERO + BSCALE*array_value           
DATE-OBS= '2016-03-18T23:22:47' /YYYY-MM-DDThh:mm:ss observation start, UT      
EXPTIME =   1.0000000000000000 /Exposure time in seconds                        EXPOSURE=   1.0000000000000000 /Exposure time in seconds                        
SET-TEMP=  -7.0000000000000000 /CCD temperature setpoint in C                   
CCD-TEMP=  -4.1500000000000004 /CCD temperature at start of exposure in C       
XPIXSZ  =   5.4000000000000004 /Pixel Width in microns (after binning)          
YPIXSZ  =   5.4000000000000004 /Pixel Height in microns (after binning)         
XBINNING=                    1 /Binning factor in width                         
YBINNING=                    1 /Binning factor in height                        
XORGSUBF=                    0 /Subframe X position in binned pixels            
YORGSUBF=                    0 /Subframe Y position in binned pixels            
READOUTM= 'Raw     ' /          Readout mode of image                           
IMAGETYP= 'LIGHT   ' /          Type of image                                   
SWCREATE= 'Celestron AstroFX V1.06' /Name of software that created the image    
COLORTYP=                    2                                                  
XBAYROFF=                    0                                                  
YBAYROFF=                    0                                                  
OBJECT  = 'test3   '                                                            
INSTRUME= 'Celestron Nightscape 8300C' /instrument or camera used               
SWOWNER = 'Dirk    ' /          Licensed owner of software                      
END                                                                                                                                                                     
*/

int NsDownload::fitsheader(fitsfile *fptr, char *fbase, struct img_params *ip)
{
    char datebuf[20] = "";
    char fbase2[64];
    struct tm t;
    gmtime_r(&ip->expdate, &t);
    strncpy(fbase2, fbase, 12);
    fbase2[12] = 0;
    snprintf(datebuf, 20, "%4d-%02d-%02dT%02d:%02d:%02d", t.tm_year + 1900, t.tm_mon + 1, t.tm_mday, t.tm_hour,
             t.tm_min, t.tm_sec);
    float xpixsize         = ip->xbinning * 4.75;
    float ypixsize         = ip->ybinning * 4.75;
    ushort zero            = 0;
    ushort two             = 2;
    int status             = 0;
    const char *image_type = NULL;
    switch (ip->frame_type) {
		case INDI::CCDChip::BIAS_FRAME:
			image_type = "BIAS";
			break;
		case INDI::CCDChip::DARK_FRAME:
			image_type = "DARK";
			break;
		case INDI::CCDChip::LIGHT_FRAME:
			image_type = "LIGHT";
			break;
		case INDI::CCDChip::FLAT_FRAME:
			image_type = "FLAT";
			break;
		default:
			DO_ERR("frame type not understood: %d", ip->frame_type);
	}

    fits_update_key(fptr, TSTRING, "DATE-OBS", datebuf, "Date of observation", &status);
    fits_update_key(fptr, TFLOAT, "EXPTIME", &ip->exp, "Exposure time in seconds", &status);
    fits_update_key(fptr, TFLOAT, "SET_TEMP", &ip->settemp, "CCD temperature setpoint in C", &status);
    fits_update_key(fptr, TFLOAT, "CCD_TEMP", &ip->acttemp, "CCD temperature at start of exposure in C", &status);
    fits_update_key(fptr, TFLOAT, "XPIXSZ", &xpixsize, "Pixel Width in microns (after binning)", &status);
    fits_update_key(fptr, TFLOAT, "YPIXSZ", &ypixsize, "Pixel Height in microns (after binning)", &status);
    fits_update_key(fptr, TSHORT, "XBINNING", &ip->xbinning, "Pixel Width in microns (after binning)", &status);
    fits_update_key(fptr, TSHORT, "YBINNING", &ip->ybinning, "Pixel Height in microns (after binning)", &status);
    fits_update_key(fptr, TUSHORT, "XORGSUBF", &zero, "Subframe X position in binned pixels", &status);
    fits_update_key(fptr, TUSHORT, "YORGSUBF", &zero, "Subframe Y position in binned pixels", &status);
    fits_update_key(fptr, TSTRING, "READOUTM", (void *)"Raw", "Readout mode of image", &status);
    fits_update_key(fptr, TSTRING, "IMAGETYP", (void *)image_type, "Type of image", &status);
    fits_update_key(fptr, TSTRING, "SWCREATE", (void *)"indi-nightscape", "Name of software that created the image", &status);
    fits_update_key(fptr, TUSHORT, "COLORTYP", &two, "Color type", &status);
    fits_update_key(fptr, TSTRING, "BAYERPAT", (void *)"BGGR", "Bayer pattern", &status);
    fits_update_key(fptr, TUSHORT, "XBAYROFF", &zero, "Bayer x offset", &status);
    fits_update_key(fptr, TUSHORT, "XBAYROFF", &zero, "Bayer y offset", &status);
    fits_update_key(fptr, TSTRING, "OBJECT", fbase, "Object being observed", &status);
    return status;
}


void NsDownload::writedownload(int pad, int cooked)
{
	char fname[64];
	char fnab[64];
	const char * extn;
	int procid = getpid();
	int nwrite = 0;
	if (cooked) {
		extn = ".fts";
	} else {
	 	extn = ".bin";	
	}
	if (retrBuf == NULL) {
		DO_DBG("%s", "no image");
		return;	
	}
	DO_INFO("done! blks %d totl %d last %d\n", retrBuf->nblks, retrBuf->nread,lastread);
	if (strlen(ctx->fbase) == 0) {
		snprintf(ctx->fbase, 64, "img_%d", procid);
	}
	if (ctx->increment) {
		snprintf(fname, 64, "%s_%d%s", ctx->fbase, ctx->imgseq, extn);
		snprintf(fnab, 64, "%s_%d", ctx->fbase, ctx->imgseq);

	} else {
		snprintf(fname, 64, "%s%s",ctx->fbase, extn);
		strcpy(fnab, ctx->fbase);
  }
	printf("%s\n",fname);

	if (pad) {
		nwrite = retrBuf->imgsz;
	} else {
		nwrite = std::min(retrBuf->nread, retrBuf->imgsz);
	}
	if (!cooked) {
		img = fopen(fname, "wb");
		if (img == NULL) {
			DO_ERR( "cannot create file %s, error %s\n", fname, strerror(errno));
		}
		
		if (img) fwrite (retrBuf->buffer, sizeof(char), nwrite , img);
		if (img) fclose(img);
	} else {
		const int max_x = get_max_x();
		int actlines = nwrite / (max_x * 2);
	  	DO_DBG("max_x: %d actlines: %d nwrite: %d imgsz: %d\n", max_x, actlines, nwrite, rd->imgsz);
		writelines = actlines;
		unsigned char * image_buf = (unsigned char *) malloc(rd->imgsz);
		const int active_x = (cn->getCamType() == kaf8300 ? KAF8300_ACTIVE_X : KAI10100_ACTIVE_X) / ctx->imgp->xbinning;
		this->copydownload(image_buf, 0, active_x, ctx->imgp->xbinning, pad, cooked);

        fitsfile *fptr;
        int status    = 0;
        long naxis    = 2, nelements;
        long naxes[2] = { active_x, actlines };
        writelines    = actlines;

        fits_create_file(&fptr, fname, &status);
        if (status != 0)
        {
            DO_ERR("attempt to create fits file '%s' failed (%d)", fname, status);
            fits_report_error(stderr, status);
            return;
        }

        fits_create_img(fptr, SHORT_IMG, naxis, naxes, &status);
        if (status != 0)
        {
            DO_ERR("attempt to create fits image (%ld, %ld, %ld) failed (%d)", naxis, naxes[0], naxes[1], status);
            fits_report_error(stderr, status);
            return;
        }

        if ((status = fitsheader(fptr, fname, ctx->imgp)) != 0)
        {
            fits_report_error(stderr, status);
            return;
        }

        nelements = naxes[0] * naxes[1];
        fits_write_img(fptr, TSHORT, 1, nelements, image_buf, &status);
        while (status != 0)
        {
			char err_buf[32];
            fits_report_error(stderr, status);
			if (ffgmsg(err_buf) == 0) { status = 0; }
            DO_ERR("attempt to write fits file %s had an error (%d) '%s'", fname, status, err_buf);
            return;
        }

        fits_close_file(fptr, &status);
        if (status != 0)
        {
            DO_ERR("attempt to close fits file %s failed (%d)", fname, status);
            fits_report_error(stderr, status);
            return;
        }

        free(image_buf);
	 DO_INFO( "wrote %d lines\n", writelines);
	}	 
}

void NsDownload::copydownload(unsigned char *buf, int xstart, int xlen, int xbin, int pad, int cooked) {
	if (cn->getCamType() == kaf8300) {
		this->copydownload_kaf8300(buf, xstart, xlen, xbin, pad, cooked);
	} else if (cn->getCamType() == kai10100) {
		this->copydownload_kai10100(buf, xstart, xlen, xbin, pad, cooked);
	} else {
		DO_ERR("bad camera type in copydownload! %d", cn->getCamType());
	}
}

void NsDownload::copydownload_kaf8300(unsigned char *buf, int xstart, int xlen, int xbin, int pad, int cooked)
{
	unsigned char linebuf[KAF8300_MAX_X * 2];
	bool forwards = true;
	bool rms = false;
	int binning = xbin;
	uint8_t * dbufp = buf;
	uint8_t * bufp;
	int nwrite = 0;
	
	if (retrBuf == NULL) {
		DO_DBG("%s", "no image");
		return;	
	}
	DO_INFO("done! blks %d totl %d last %d\n", retrBuf->nblks, retrBuf->nread,lastread);
			
	if (!cooked) {
		if (pad) {
			nwrite = retrBuf->imgsz;
		} else {
			nwrite = retrBuf->nread;
		}
		memcpy (dbufp, retrBuf->buffer, nwrite);
	} else {
		//int actlines = nwrite / (KAF8300_MAX_X*2);
	  //int rem = nwrite % (KAF8300_MAX_X*2);
	  nwrite = retrBuf->nread;
		int nwriteleft = nwrite;
		if (forwards) {
			bufp = retrBuf->buffer;
			//dbufp =  (dbufp +(KAF8300_ACTIVE_X*2*IMG_Y)) - (KAF8300_ACTIVE_X*2);
		} else {
			bufp = (retrBuf->buffer + nwrite) - (KAF8300_MAX_X*2);
			dbufp =  (dbufp +(KAF8300_ACTIVE_X*2*KAF8300_IMG_MAX_Y)) - (KAF8300_ACTIVE_X*2);
		}
		writelines = 0;
	  while (nwriteleft >= (KAF8300_MAX_X*2)) {
	  //swab (bufp + (KAF8300_POSTAMBLE*2),  linebuf, KAF8300_ACTIVE_X*2 );
	   //memcpy (dbufp, linebuf, KAF8300_ACTIVE_X*2);
	    if (binning > 1) {
	    	uint8_t * lbufp = bufp + (KAF8300_POSTAMBLE*2) + xstart*2;
	    	int len = xlen * 2;
	    	int linelen = 0;
	    	while (len > 0) {
	    		short px[4];
	    		long long pxsq =0;
	    		long pxav =0;
	    		short pxa;
	    		memcpy(px, lbufp,binning * 2);
	    		for (int a = 0; a < binning; a++) {
	    			pxav += px[a];
	    			pxsq	+= px[a]*px[a];
	    		}
	    		pxav /= binning;
	    		pxsq /= binning;
	    		if(rms) {
	    			pxa = round(sqrt((double)	pxsq));
	    			
	    		} else {
	    			pxa = pxav;	
	    		}
	    		memcpy (linebuf + linelen, &pxa, 2);
	    		linelen += 2;
	    		lbufp += 2*binning;
	    		len -= 2* binning;
				}
				memcpy (dbufp, linebuf, (xlen*2)/binning);

	    } else {
	  		memcpy (dbufp, bufp + (KAF8300_POSTAMBLE*2) + xstart*2, xlen * 2 ); //KAF8300_ACTIVE_X*2 );
	  	}
			if (forwards) { 
				bufp +=  KAF8300_MAX_X*2;
				dbufp +=(xlen*2)/binning;
			} else {
				bufp -=  KAF8300_MAX_X*2;
				dbufp-=(xlen*2)/binning;
			}
			nwriteleft -= KAF8300_MAX_X*2;
			writelines++;
	  }
	 DO_INFO( "wrote %d lines\n", writelines);
	}	 
}

void NsDownload::copydownload_kai10100(unsigned char *image, int xstart, int xlen, int xbin, int pad, int cooked) {
	uint8_t * img_dest = image;
	uint8_t * img_src;
	int nwrite = 0;
	
	if (retrBuf == NULL) {
		DO_DBG("%s", "no image");
		return;	
	}
	DO_DBG("copydownload image: %p xstart: %d xlen: %d xbin: %d pad: %d cooked: %d\n", image, xstart, xlen, xbin, pad, cooked);
	DO_INFO("done! blks %d totl %d last %d\n", retrBuf->nblks, retrBuf->nread,lastread);
			
	if (!cooked) {
		if (pad) {
			nwrite = retrBuf->imgsz;
		} else {
			nwrite = retrBuf->nread;
		}
		memcpy (img_dest, retrBuf->buffer, nwrite);
	} else {
		const int max_x = get_max_x();
		const int active_y = xbin == 1 ? KAI10100_ACTIVE_Y : xbin == 2 ? KAI10100_HALF_Y : KAI10100_QUARTER_Y;

		int actual_lines = retrBuf->nread / (max_x * 2);
		DO_DBG("actual_lines: %d nread: %d max_x: %d\n", actual_lines, retrBuf->nread, max_x);
		if (actual_lines > active_y) {actual_lines = active_y;}

	  	nwrite = retrBuf->nread;
		writelines = 0;
	    if (xbin == 1) {
			// At binning 1, the image is downloaded in 4 groups which need to be interlaced
			// together to form the final image.
			for (int offset = 0; offset < 4; ++offset) {
				img_src = retrBuf->buffer + (offset * KAI10100_INTERLACE * max_x * 2) + ((KAI10100_X_PREAMBLE + xstart) * 2);
				img_dest = image + (offset * xlen * 2);
				for (int line = 0; line < (actual_lines / 4); ++line) {
					memcpy(img_dest, img_src, xlen * 2);
					img_src += max_x * 2;
					img_dest += xlen * 2 * 4;
					writelines++;
				}
			}
	    } else {
			img_src = retrBuf->buffer + ((KAI10100_X_PREAMBLE / xbin + xstart) * 2);
			img_dest = image;
			for (int line = 0; line < actual_lines; ++line) {
				memcpy(img_dest, img_src, xlen * 2);
				img_src +=  max_x * 2;
				img_dest += (xlen * 2);
				writelines++;
			}
	  	
		  }
		DO_INFO( "wrote %d lines\n", writelines);
	}	 
}


int NsDownload::purgedownload() 
{
		int rc2;
		rc2 = cn->readData(rd->buffer, rd->bufsiz);
		if (rc2 < 0 ) {
			DO_ERR( "purge: unable to read: %d \n", rc2);
			return (-1);
		}
		if (rc2 > 0) {
			DO_ERR("purge: spare read %d\n", rc2);
		  rc2 = cn->purgeData();
		  if (rc2 < 0 ) {
		    DO_ERR( "unable to purge: %d\n", rc2);
				return (-1);
			}
		}	
		return 0;
}


int NsDownload::fulldownload() 
{
		int rc2;
		
		rc2 = cn->readData(rd->buffer+rd->nread, rd->bufsiz - rd->nread);
		if (rc2 < 0 ) {
			DO_ERR( "unable to read: %d\n", rc2);
			return (-1);
		}
		if (rc2 > 0) {
			DO_INFO("read %d\n", rc2);
		  rd->nread += rc2;
		  rd->nblks += rc2/65536;
		  DO_INFO("read %d tot %d\n", rc2, rd->nread);

		}	
		return rc2;
}


void NsDownload::initdownload()
{
	long imgszmax;
	if (cn->getCamType() == kaf8300) {
		imgszmax = KAF8300_MAX_X*KAF8300_ACTIVE_Y*2 + DEFAULT_CHUNK_SIZE;
	} else {
		imgszmax = KAI10100_MAX_X*KAI10100_MAX_Y*2 + DEFAULT_CHUNK_SIZE;
	}
		readdone = 0;
		rd->nread = 0;
		if(!rd->buffer) {
			rd->buffer = (unsigned char *)malloc(imgszmax);
		}
		memset(rd->buffer, 0, imgszmax);

		rd->bufsiz = imgszmax;
		rd->nblks = 0;	
}


// static void  download_thread(NsDownload * d) {
//	d->trun();	
//}


void NsDownload::trun()
{

	int maxxfer;
  int zeroes = 0;
	//struct ftdi_context * ftdid = cn->getDataChannel();
	maxxfer = cn->getMaxXfer();

	do  {
		DO_INFO("%s\n", "initdownload");
		std::unique_lock<std::mutex> ulock(mutx);

		while(!do_download && !interrupted) go_download.wait(ulock);
		if(interrupted) break;
		initdownload();

	
		if (do_download && !in_download) {
			//ftdi_usb_close(ftdid);
			//maxxfer = opendownload(ftdid, ctx->dev);
			in_download = 1;
			ctx->imgseq++;
			zeroes = 0;
		}
	  while (in_download && !interrupted) {
	  	//int rc2= cn->setDataRts();;
    	//if (rc2  < 0) {
      //  		DO_ERR( "unable to set rts: %d\n", rc2);
    	//}	
    	int down = 0;
	    if (zero_reads > 1) 
	  		down = fulldownload();
	  	else
	  		down = downloader();
	  	if (down < 0) {
	  		DO_ERR( "unable to read download: %d\n", down);
	  		do_download = 0;
	  		in_download = 0;
	  		continue;
	  	}
	  	if (rd->nread < rd->imgsz) {
	  		if (down == 0 && rd->nread > 0) {
    			zeroes++;
	  		}
	  		if (zeroes < zero_reads) continue;
	  	}
	  	int pad = 0;
	  	if (rd->nread != rd->imgsz) {
			const int max_x = get_max_x();
	  		int actlines = rd->nread / (max_x*2);
	  		int rem = rd->nread % (max_x*2);
	  		DO_INFO( "siz %d read %d act lines %d rem %d\n",  rd->imgsz, rd->nread, actlines, rem);
	  		if (rd->imgsz - rd->nread < max_x * 5) {
	  			pad = 1;
	  		}
	    }	
	    lastread = down;
	    	   // IDLog("foop\n");

	    if (zero_reads > 1) {
	    	rb = rdd;
	    	retrBuf = &rb;
	    	rd->buffer = NULL;
	    }
	    //IDLog("retr %p buf %p \n", retrBuf, rb.buffer);
	    if(write_it) writedownload(pad, 1);
			
	  	do_download = 0;
	  	in_download = 0;
	  }
	  if (!in_download && !do_download) {
			initdownload();  	
	  	purgedownload ();
	  }
  } 	while(ctx->nexp >= ctx->imgseq && !interrupted);
  		DO_DBG("%s\n", "thread done");

}


void NsDownload::startThread(void) {
			interrupted = 0;
			sched_param sch_params;
      sch_params.sched_priority = 3;
    	downthread = new std::thread (&NsDownload::trun, this); //(&NsDownload::trun, this);
			pthread_setschedparam(downthread->native_handle(), SCHED_FIFO, &sch_params);		

}

void NsDownload::stopThread(void) {
	setInterrupted();
	downthread->join();
}
