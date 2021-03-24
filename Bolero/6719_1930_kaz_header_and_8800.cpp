#include <pragmas_2411.h>



unsigned char *data = (unsigned char*)malloc(2 * 44100);
signed short *sound = (signed short*)malloc(2 * 44100);

void main(void){
	FILE *fp;
	int i;
	int data_no;
	fopen_s(&fp, "6616_London_g-1.wav", "rb");
	int totals1,totals2;
	if (fp != NULL) {

		fread(data, 1, 44, fp);
//---------------------------------------------------------------------------
		totals1 = data[4] + 256 * (data[5] + 256 * (data[6] + 256 * data[7]));
		int fmt_ch = data[16] + 256 * (data[17] + 256 * (data[18] + 256 * data[19]));
		int s_rate = data[24] + 256 * (data[25] + 256 * (data[26] + 256 * data[27]));
		int speed = data[28] + 256 * (data[29] + 256 * (data[30] + 256 * data[31]));
		int block_size = data[32] + 256 * data[33];
		int bit_no = data[34] + 256 * data[35];
		totals2 = data[40] + 256 * (data[41] + 256 * (data[42] + 256 * data[43]));

		printf_s("\n  0- 3 RIFF                         %c%c%c%c", data[0], data[1], data[2], data[3]);
		printf_s("\n  4- 7 file size                    %5x %5d", totals1, totals1);
		printf_s("\n  8-12 WAV                          %c%c%c%c", data[8], data[9], data[10], data[11]);
		printf_s("\n 12-15 fmt                          %c%c%c%c", data[12], data[13], data[14], data[15]);
		printf_s("\n 16-19 size of fmt chunk            %5d", fmt_ch);
		printf_s("\n 20-21 format ID                    %5d", data[20]+256*data[21]);
		printf_s("\n 22-23 channel number               %5d", data[22]+256*data[23]);
		printf_s("\n 24-27 sampling rate(Hz)            %5x  %5d", s_rate,s_rate);
		printf_s("\n 28-31 speed(= byte/s)              %5x  %5d", speed,speed);
		printf_s("\n 32-33 block size(=ch*bit/8byte)    %5d", block_size);
		printf_s("\n 34-35 bit number(8,16,24,32bit)    %5d",bit_no);
		printf_s("\n 36-39 data                         %c%c%c%c", data[36], data[37], data[38], data[39]);
		printf_s("\n 40-43 data size(byte)              %5x  %5d", totals2, totals2);

 //-----------------------------------------------------------------------------------
		printf_s("\n\n read actual data from here\n             ");
		IplImage *grafs;
		cvNamedWindow("grafs",1);
		cvMoveWindow("grafs", 400,0);

		grafs = cvCreateImage(cvSize(600, 600), IPL_DEPTH_8U, 3);


		for (int loop = 0; loop < 60; loop++){
			cvRectangle(grafs, cvPoint(0, 0), cvPoint(600, 600), CV_RGB(255, 255, 255), -1, 8, 0);
			cvRectangle(grafs, cvPoint(50, 50), cvPoint(550, 250), CV_RGB(0, 0, 0), 2, 8, 0);
			cvLine(grafs, cvPoint(50, 150), cvPoint(550, 150), CV_RGB(0, 0, 0), 2, 8, 0);

			cvRectangle(grafs, cvPoint(50, 300), cvPoint(550, 500), CV_RGB(0, 0, 0), 2, 8, 0);
			cvLine(grafs, cvPoint(50, 400), cvPoint(550, 400), CV_RGB(0, 0, 0), 2, 8, 0);
			
			data_no = fread(sound, 2, 2 * 10000 * 2, fp);
			if (data_no < 2 * 8820 * 2) {
				break;
			}
			for (int xx = 0; xx < data_no / 2 - 2; xx++) {
				cvLine(grafs, cvPoint(50 + xx / 40, 150 + sound[2 * xx]/10),
					cvPoint(50 + (xx + 1) / 40, 150 + sound[2 * xx + 2]/10), CV_RGB(255, 0, 0), 1, 8, 0);
				cvLine(grafs, cvPoint(50 + xx / 40,400 + sound[2 * xx + 1]/10),
					cvPoint(50 + (xx + 1) / 40, 400 + sound[2 * xx + 3]/10), CV_RGB(0, 255, 0), 1, 8, 0);
			}
			cvShowImage("grafs", grafs);
			cvWaitKey(50);
		}
		}
	cvWaitKey(0);
		fclose(fp);
	}
