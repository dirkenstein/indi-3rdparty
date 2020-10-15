#ifndef __NS_CHANNEL_H__
#define __NS_CHANNEL_H__
#include <stdlib.h>
#include <libftdi1/ftdi.h>
#define DEFAULT_OLD_CHUNK_SIZE  63448
#define DEFAULT_CHUNK_SIZE 65536

typedef enum {
	none = 0,
	kaf8300,
	kai10100
} camera_t;

class NsChannel {
	public:
		NsChannel() {
			maxxfer = 0;
			opened = 0;
			camnum = 0;
			cam_type = none;
		}
		NsChannel(int cam) {
			camnum = cam;
			maxxfer = 0;
			opened = 0;
			cam_type = none;
		}
		
		virtual ~NsChannel() { if (opened) close(); };
		int open();
		int getMaxXfer();
		camera_t getCamType();
		void setCamType(camera_t camType);

	  virtual int close();
		virtual int readCommand(unsigned char * buf, size_t n) = 0;
		virtual int writeCommand(const unsigned char * buf, size_t n) = 0;
		virtual int readData(unsigned char * buf, size_t n)= 0;
		virtual int purgeData(void)= 0;
		virtual int setDataRts(void)= 0;
		virtual int resetcontrol (void)= 0;

	protected:
		virtual int opencontrol (void)= 0;

		virtual int opendownload(void)= 0;
		virtual int scan(void)= 0;
		static const int vid = 0x19b4;
    static const int pid = 0x0065;
		unsigned camnum;
		int maxxfer;
		unsigned ndevs;
		bool opened;
		int thedev;
	camera_t cam_type;
	

};

#endif
