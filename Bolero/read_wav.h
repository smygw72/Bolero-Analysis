#include <string>
#include <string.h>
#include <iostream>
using namespace std;

namespace sl
{
	// struct -----------------------------------------------------
	struct Wav // basic WAVE file structure
	{
	public:
		int fs; // Samples per Second
		int bits; //  Bits per Sample
		int channel; // Sound Data Channels
		int length; // Sound Data Length
		vector<vector<double>> data;
	};
	
	// function --------------------------------------------------
	Wav wavread(string wavPath);
	void wavwrite(string wavPath, Wav wavData);
	Wav stereo2mono(Wav wavData);
	Wav mono2stereo(Wav wavData, double delayTime = 0.005);
	void showWavData(Wav wavData);
	double sinc(double);
	
	Wav wavread(string wavPath)
	{
		/*MATLAB like wavread function*/
		// input : WAVE file absolute path
		// output : Wav struct data defined above
		const char *wavPath_c = wavPath.c_str();
		Wav output;

		FILE *fp;
		errno_t error;
		int n, m;

		char riff_chunkID[4];
		long riff_chunkSize;
		char riff_formType[4];
		char fmt_chunkID[4];
		long fmt_chunkSize;
		short fmt_waveFormatType;
		short fmt_channel;
		long fmt_samplesPerSec;
		long fmt_bytesPerSec;
		short fmt_blockSize;
		short fmt_bitsPerSample;
		char data_chunkID[4];
		long data_chunkSize;
		short data;

		long fmt_extend;
		char find_data_char[2]; // To find "data"
		char find_data_char_old[2]; // To find "data"
		int dataBlockSize; // data size per sample

		if (error = fopen_s(&fp, wavPath_c, "rb") != 0){
			cout << "ERROR : Cannot read wave file. Check file path." << endl;
			output.bits = -1;
			output.channel = -1;
			output.data.resize(1);
			output.data[0].resize(1);
			output.fs = -1;
			output.length = -1;
			return output;
		}

		// RIFF
		fread(riff_chunkID, 1, 4, fp);
		fread(&riff_chunkSize, 4, 1, fp);
		fread(riff_formType, 1, 4, fp);

		// fmt chunk
		fread(fmt_chunkID, 1, 4, fp);
		fread(&fmt_chunkSize, 4, 1, fp); // Linear PCM = 16
		fread(&fmt_waveFormatType, 2, 1, fp);
		fread(&fmt_channel, 2, 1, fp);
		fread(&fmt_samplesPerSec, 4, 1, fp);
		fread(&fmt_bytesPerSec, 4, 1, fp);
		fread(&fmt_blockSize, 2, 1, fp);
		fread(&fmt_bitsPerSample, 2, 1, fp);
		if (fmt_chunkSize != 16){ // if WAVE file is not linear PCM
			cout << "WARNING : This file is not linear PCM" << endl;
			cout << "There is a possibility to be malfunctioning." << endl;
			fread(&fmt_extend, fmt_chunkSize -16, 1, fp);
		}

		// data chunk
		fread(data_chunkID, 1, 4, fp);
		if (memcmp(data_chunkID,"data", sizeof(data_chunkID))!=0){// find "data" in the binary data
			fread(find_data_char_old, 1, 2, fp);
			while (1){
				fread(find_data_char, 1, 2, fp);
				if (memcmp(find_data_char_old,"da", sizeof(find_data_char_old))==0 && memcmp(find_data_char,"ta", sizeof(find_data_char))==0) // find "data"
					break;
				find_data_char_old[0] = find_data_char[0];
				find_data_char_old[1] = find_data_char[1];
			}
		}
		fread(&data_chunkSize, 4, 1, fp);

		output.fs = fmt_samplesPerSec;
		output.bits = fmt_bitsPerSample;
		output.length = data_chunkSize / (fmt_bitsPerSample / 8) / fmt_channel;
		output.channel = fmt_channel;
		output.data.resize(output.length);
		for (int i = 0; i < output.data.size(); ++i) {
			output.data[i].resize(output.channel);
		}

		// read wave data
		dataBlockSize = fmt_blockSize / fmt_channel;
		if (fmt_channel == 1){
			for (n = 0; n < output.length; n++){
				fread(&data, dataBlockSize, 1, fp);
				output.data[n][0] = (double) data / 32768.0; // normalize in range from -1 to 1
			}
		}else if(fmt_channel == 2){
			for (n = 0; n < output.length; n++){
				for (m=0; m<2; m++){
					fread(&data, dataBlockSize, 1, fp);
					output.data[n][m] = (double) data / 32768.0; // normalize in range from -1 to 1
				}
			}
		}

		fclose(fp);
		return output;
	}

	void wavwrite(string wavPath, Wav wavData)
	{
		/*Write WAVE file with linear PCM, 16bit*/
		// argument 1 : wave file path
		// argument 2 : wave file data
		const char *wavPath_c = wavPath.c_str();
		FILE *fp;
		errno_t error;
		int n, m;

		char riff_chunkID[4];
		long riff_chunkSize;
		char riff_formType[4];
		char fmt_chunkID[4];
		long fmt_chunkSize;
		short fmt_waveFormatType;
		short fmt_channel;
		long fmt_samplesPerSec;
		long fmt_bytesPerSec;
		short fmt_blockSize;
		short fmt_bitsPerSample;
		char data_chunkID[4];
		long data_chunkSize;
		short data;
		vector<vector<double>> mat;
		mat.resize(wavData.data.size());
		for (int i = 0; i < mat.size(); ++i) {
			for (int j = 0; j < wavData.data[i].size();++j) {
				mat[i].push_back(1);
			}
		}
		int dataBlockSize;

		riff_chunkID[0] = 'R';
		riff_chunkID[1] = 'I';
		riff_chunkID[2] = 'F';
		riff_chunkID[3] = 'F';
		riff_chunkSize = 36 + wavData.length * wavData.bits/8 * wavData.channel;
		riff_formType[0] = 'W';
		riff_formType[1] = 'A';
		riff_formType[2] = 'V';
		riff_formType[3] = 'E';

		fmt_chunkID[0] = 'f';
		fmt_chunkID[1] = 'm';
		fmt_chunkID[2] = 't';
		fmt_chunkID[3] = ' ';
		fmt_chunkSize = 16; // linear PCA : 16
		fmt_waveFormatType = 1;
		fmt_channel = wavData.channel;
		fmt_samplesPerSec = wavData.fs;
		fmt_bytesPerSec = wavData.fs * wavData.bits/8 * wavData.channel;
		fmt_blockSize = wavData.bits/8 * wavData.channel;
		fmt_bitsPerSample = wavData.bits;

		data_chunkID[0] = 'd';
		data_chunkID[1] = 'a';
		data_chunkID[2] = 't';
		data_chunkID[3] = 'a';
		data_chunkSize = wavData.length * wavData.bits/8 * wavData.channel;

		if (error = fopen_s(&fp, wavPath_c, "wb") != 0){
			cout << "ERROR : Cannot open " << wavPath << endl;
			return;
		}

		// RIFF
		fwrite(riff_chunkID, 1, 4, fp);
		fwrite(&riff_chunkSize, 4, 1, fp);
		fwrite(riff_formType, 1, 4, fp);

		// fmt chunk
		fwrite(fmt_chunkID, 1, 4, fp);
		fwrite(&fmt_chunkSize, 4, 1, fp); // Linear PCM = 16
		fwrite(&fmt_waveFormatType, 2, 1, fp);
		fwrite(&fmt_channel, 2, 1, fp);
		fwrite(&fmt_samplesPerSec, 4, 1, fp);
		fwrite(&fmt_bytesPerSec, 4, 1, fp);
		fwrite(&fmt_blockSize, 2, 1, fp);
		fwrite(&fmt_bitsPerSample, 2, 1, fp);

		// data chunk
		fwrite(data_chunkID, 1, 4, fp);
		fwrite(&data_chunkSize, 4, 1, fp);

		dataBlockSize = fmt_blockSize / fmt_channel;
		for (int i = 0; i < mat.size(); ++i) {
			for (int j = 0; j < mat[i].size(); ++j) {
				mat[i][j] = (mat[i][j] + wavData.data[i][j]) / 2.0 * 65536;
			}
		}
		if (fmt_channel == 1){
			for (n=0; n< wavData.length; n++){
				if (mat[n][0]> 65535.0){
					mat[n][0] = 65535.0; // clipping
				}else if (mat[n][0] < 0){
					mat[n][0] = 0.0; // clipping
				}
				data = (short) (mat[n][0] + 0.5) -32768; // round off and adjust off-set
				fwrite(&data, dataBlockSize, 1, fp);
			}
		}else if(fmt_channel == 2){
			for (n=0; n< wavData.length; n++){
				for (m=0; m<2; m++){
					if (mat[n][m] > 65535.0){
						mat[n][m] = 65535.0; // clipping
					}else if (mat[n][m] < 0){
						mat[n][m] = 0.0; // clipping
					}
					data = (short) (mat[n][m] + 0.5) -32768; // round off and adjust off-set
					fwrite(&data, dataBlockSize, 1, fp);
				}
			}
		}

		fclose(fp);
	}

	Wav stereo2mono(Wav wavData)
	{
		/*Convert WAVE data from stereo into monaural signal*/
		Wav output;
		if (wavData.channel != 2){
			cout << "ERROR : The number of channels must be 2." << endl;
			return wavData;
		}
		output.bits = wavData.bits;
		output.channel = 1;
		output.fs = wavData.fs;
		output.length = wavData.length;
		for (int i = 0; i < output.data.size(); ++i) {
			output.data[i][0] = (wavData.data[i][0] + wavData.data[i][1]) / 2.0;
			output.data[i].resize(0);
		}
		return output;
	}

	Wav mono2stereo(Wav wavData, double delayTime)
	{
		/*Convert WAVE data from monaural signal into stereo*/
		// quasi-stereo

		Wav output;
		vector<vector<double>> mat;
		int tempTime;
		if (wavData.channel != 1){
			cout << "ERROR : The number of channels must be 1." << endl;
			return output;
		}
		output.bits = wavData.bits;
		output.channel = 2;
		output.fs = wavData.fs;
		output.length = wavData.length;
		mat.resize(wavData.length);
		for (int i = 0; i < wavData.length; ++i) {
			for (int j = 0; j < 2; ++j) mat[i].push_back(0);
			mat[i][0] = wavData.data[i][0];
			mat[i][1] = wavData.data[i][0];
		}

		delayTime *= wavData.fs;

		for (int i=0; i<wavData.length; i++){
			tempTime = i - (int) delayTime;
			if (tempTime >= 0){
				mat[i][0] += wavData.data[tempTime][0];
				mat[i][0] -= wavData.data[tempTime][0];
			}
		}
		output.data = mat;
		return output;
	}

	void showWavData(Wav wavData)
	{
		/*Print WAVE file data*/
		cout << "Bits : " << wavData.bits << endl;
		cout << "Channel : " << wavData.channel << endl;
		cout << "Sampling rate : " << wavData.fs << endl;
		cout << "Data length : " << wavData.length << endl;
	}

	int log2(int num)
	{
		if (num<0){
			return -1;
		}
		return log(num) / log(2);
	}

	double sinc(double x)
	{
		double y;
		if(x == 0.0){
			y = 1.0;
		}else{
			y = sin(x) / x;
		}
		return y;
	}

}