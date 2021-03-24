#include <iostream>
#include <fstream>
#include <Eigen/Core>
#include <opencv2/opencv.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <complex.h>
#define _USE_MATH_DEFINES 
#include <math.h>
#include <fftw3.h>
#include "read_wav.h"

#define BOLERO 1 // ボレロ
#define SIN_WAVE 0 // 校正用
#define CF 1 // 相関関数の計算の有無

// Global変数
sl::Wav wav;
uint waveSize;
double timesX;

// 縦軸のスケール修正
const double SCALE = 4.5315*pow(10, -5);
// ただし、この修正では対数値がマイナスになってしまうためデータ値全体を10^5倍する

//マウス入力用のパラメータ
struct mouseParam {
	int x;
	int y;
	int event;
	int flags;
};

//コールバック関数
void CallBackMouseFunc(int eventType, int x, int y, int flags, void* userdata)
{
	mouseParam *ptr = static_cast<mouseParam*> (userdata);

	ptr->x = x;
	ptr->y = y;
	ptr->event = eventType;
	ptr->flags = flags;
}

/* 相関関数の描画 */
void drawGraph(vector<double> _data,const string str) 
{
	cv::Mat img(600, 1200, CV_8UC3, cv::Scalar(255, 255, 255));
	const int startX = 50, startY = 50, endX = 1150, endY = 550;

	cv::line(img, cv::Point(startX, startY), cv::Point(startX, endY), cv::Scalar(0, 0, 0), 4, 4);
	cv::line(img, cv::Point(startX, startY), cv::Point(endX, startY), cv::Scalar(0, 0, 0), 4, 4);

	double timesX = (double)(endX - startX) / (log10(_data.size()));
	double timesY = (double)(endY - startY) / 1.0;

	for (int i = 1; i < _data.size(); ++i) {
		if (_data[i] >= 0) {

			double noramlizedDataX = startX + log10(i) * timesX;
			double normalizedDataY = startY + _data[i] * timesY;

			cv::line(img, cv::Point(noramlizedDataX, startY), cv::Point(noramlizedDataX, normalizedDataY), cv::Scalar(0, 0, 0), 1, 4);
		}
	}

	cv::flip(img, img, 0);

	// x
	cv::putText(img, "100", cv::Point(startX + 2 * timesX, 570), cv::FONT_HERSHEY_TRIPLEX, 0.5, cv::Scalar(0, 0, 0), 1, CV_AA);
	cv::putText(img, "500", cv::Point(startX + log10(500) * timesX, 570), cv::FONT_HERSHEY_TRIPLEX, 0.5, cv::Scalar(0, 0, 0), 1, CV_AA);
	cv::putText(img, "1000", cv::Point(startX + 3 * timesX, 570), cv::FONT_HERSHEY_TRIPLEX, 0.5, cv::Scalar(0, 0, 0), 1, CV_AA);
	cv::putText(img, "10000", cv::Point(startX + 4 * timesX, 570), cv::FONT_HERSHEY_TRIPLEX, 0.5, cv::Scalar(0, 0, 0), 1, CV_AA);

	cv::imshow(str, img);
}

/* 自己相関関数 */
void calculateAcf(const vector<vector<double>> _input) 
{
	/* フーリエ変換 */
	vector<double> dataY;

	// FFTW
	fftwf_complex *in, *out;
	fftwf_plan plan;

	// メモリ確保
	in = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*_input.size());
	out = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*_input.size());

	// 基本角周波数
	double w = (double)wav.fs / (double)_input.size();

	// 初期化
	for (int i = 0; i < _input.size(); ++i)
	{
		in[i][0] = _input[i][0];
		in[i][1] = 0.0;
		out[i][0] = 0.0;
		out[i][1] = 0.0;
	}

	plan = fftwf_plan_dft_1d(_input.size(), in, out, FFTW_FORWARD, FFTW_ESTIMATE);

	fftwf_execute(plan);
	fftwf_destroy_plan(plan);

	for (int i = 0; i < (_input.size() / 2); i++)
	{
		double re = out[i][0];
		double im = out[i][1];
		double mag = sqrt(re*re + im*im);
		dataY.push_back(mag);
	}

	// 解放
	fftwf_free(in);
	fftwf_free(out);

	fftwf_complex *in_inverse, *out_inverse;
	fftwf_plan plan_inverse;

	in_inverse = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*dataY.size());
	out_inverse = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*dataY.size());

	for (int i = 0; i < dataY.size(); ++i) {
		in_inverse[i][0] = dataY[i];
		in_inverse[i][1] = 0.0;
		out_inverse[i][0] = 0.0;
		out_inverse[i][1] = 0.0;
	}

	plan_inverse = fftwf_plan_dft_1d((dataY.size() / 2), in_inverse, out_inverse, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftwf_execute(plan_inverse);
	cout << "逆フーリエ変換実行" << endl;

	vector<double> acfY;
	double maxY = 0.0;
	for (int i = 0; i < (_input.size() / 2); i++)
	{
		double re = out_inverse[i][0];         // 複素数の実数部
		double im = out_inverse[i][1];         // 複素数の虚数部
		//cout << im << endl;
		double mag = sqrt(re*re + im*im);      // 大きさを計算
		acfY.push_back(mag);
		if (maxY < mag) maxY = mag;
	}

	for (int i = 0; i < acfY.size(); ++i) {
		acfY[i] /= maxY;
		//cout << acfY[i] << endl;
	}
	fftwf_destroy_plan(plan_inverse);
	drawGraph(acfY, "自己相関関数");
}

/* 相互相関関数 */
void calculateCcf(vector<vector<double>> _input)
{
	/* フーリエ変換 */
	vector<double> data1Re, data1Im, data2Re, data2Im;

	// FFTW
	fftwf_complex *in1, *out1, *in2, *out2;
	fftwf_plan plan1, plan2;

	// メモリ確保
	in1 = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*_input.size());
	out1 = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*_input.size());
	in2 = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*_input.size());
	out2 = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*_input.size());

	// 基本角周波数
	double w = (double)wav.fs / (double)_input.size();

	// 初期化
	for (int i = 0; i < _input.size(); ++i)
	{
		in1[i][0] = _input[i][0];
		in1[i][1] = 0.0;
		out1[i][0] = 0.0;
		out1[i][1] = 0.0;
		in2[i][0] = _input[i][1];
		in2[i][1] = 0.0;
		out2[i][0] = 0.0;
		out2[i][1] = 0.0;
	}

	plan1 = fftwf_plan_dft_1d(_input.size(), in1, out1, FFTW_FORWARD, FFTW_ESTIMATE);
	plan2 = fftwf_plan_dft_1d(_input.size(), in2, out2, FFTW_FORWARD, FFTW_ESTIMATE);

	fftwf_execute(plan1);
	fftwf_destroy_plan(plan1);
	fftwf_execute(plan2);
	fftwf_destroy_plan(plan2);

	for (int i = 0; i < (_input.size() / 2); i++)
	{
		double re1 = out1[i][0];
		double im1 = out1[i][1];
		double re2 = out2[i][0];
		double im2 = out2[i][1];

		data1Re.push_back(re1);
		data1Im.push_back(im1);
		data2Re.push_back(re2);
		data2Im.push_back(im2);
	}

	// 解放
	fftwf_free(in1);
	fftwf_free(out1);
	fftwf_free(in2);
	fftwf_free(out2);

	/* 逆フーリエ変換 */
	fftwf_complex *in_inverse, *out_inverse;
	fftwf_plan plan_inverse;

	// メモリ確保
	in_inverse = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*data1Re.size());
	out_inverse = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*data1Re.size());

	// 初期化
	for (int i = 0; i < data1Re.size(); ++i)
	{
		in_inverse[i][0] = data1Re[i] * data2Re[i] + data1Im[i] * data2Im[i];
		in_inverse[i][1] = data1Im[i] * data2Re[i] - data1Re[i] * data2Im[i];
		out_inverse[i][0] = 0.0;
		out_inverse[i][1] = 0.0;
	}

	plan_inverse = fftwf_plan_dft_1d(data1Re.size(), in_inverse, out_inverse, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftwf_execute(plan_inverse);
	fftwf_destroy_plan(plan_inverse);

	vector<double> ccf;
	double maxY = 0.0;
	int index = -1;
	for (int i = 0; i < (data1Re.size() / 2); i++)
	{
		double re = out_inverse[i][0];
		double im = out_inverse[i][1];
		double mag = sqrt(re*re + im*im);
		ccf.push_back(mag);
		if (mag > maxY) {
			maxY = mag;
			index = i;
		}
	}
	cout << index << endl;
	fftwf_free(in_inverse);
	fftwf_free(out_inverse);

	for (int i = 0; i < ccf.size(); ++i) {
		ccf[i] /= maxY;

	}
	drawGraph(ccf, "相互相関関数");
}

/* 周波数応答 */
void DrawFrequency(cv::Mat& _img, const vector<double> _dataX, const vector<double> _dataY) 
{
	_img = cv::Mat(600, 1200, CV_8UC3, cv::Scalar(255, 255, 255));

	// グラフ描画範囲 原点(50,50)
	const int startX = 50, startY = 50, endX = 1150, endY = 550;

	// 軸
	cv::line(_img, cv::Point(startX, startY), cv::Point(startX, endY), cv::Scalar(0, 0, 0), 4, 4);
	cv::line(_img, cv::Point(startX, startY), cv::Point(endX, startY), cv::Scalar(0, 0, 0), 4, 4);

	// 1番大きい値を探す
	int index = 0;
	double maxValueY = 0.0;
	int maxValueX = 0;

	for (int i = 0; i < _dataY.size(); ++i) {
		if (maxValueY < _dataY[i]) {
			maxValueY = _dataY[i];
		}
	}
	// 正規化
	if (_dataX.size()!=0 && _dataX.back() > 2) {
		timesX = (double)(endX - startX) / (_dataX.back() - 2);
	}
	double timesY = (double)(endY - startY) / maxValueY;

	// 表示
	for (int i = 0; i < _dataY.size(); ++i) {

		if (_dataX[i] >= 2 && _dataY[i] >= 0) {
		
			// 低周波数領域は興味ないので100Hz以上を表示するようにする(1/100倍(-2)する)
			double noramlizedDataX = startX + (_dataX[i] - 2) * timesX;
			double normalizedDataY = startY + _dataY[i]* timesY;
			
			if (_dataY[i] == maxValueY) {
				//cv::line(frequency, cv::Point(noramlizedDataX, startY), cv::Point(noramlizedDataX, normalizedDataY), cv::Scalar(0, 0, 255), 4, 4);
				//cv::line(frequency, cv::Point(50, normalizedDataY), cv::Point(1150, normalizedDataY), cv::Scalar(0, 0, 255), 4, 4);

				maxValueX = pow(10, _dataX[i]);
				cv::flip(_img, _img, 0);

				cv::putText(_img, to_string(maxValueX), cv::Point(noramlizedDataX - 10, 570), cv::FONT_HERSHEY_TRIPLEX, 0.5, cv::Scalar(0, 0, 255), 1, CV_AA);
				cv::putText(_img, to_string(waveSize), cv::Point(30, 550-normalizedDataY), cv::FONT_HERSHEY_TRIPLEX, 0.5, cv::Scalar(0, 0, 255), 1, CV_AA);

				cv::flip(_img, _img, 0);
			}
			else {
				cv::line(_img, cv::Point(noramlizedDataX, startY), cv::Point(noramlizedDataX, normalizedDataY), cv::Scalar(0, 0, 0), 1, 4);
			}
		}
	}

	// 反転
	cv::flip(_img, _img, 0);

	/* 目盛 */
	cv::putText(_img, "[Hz]", cv::Point(1150, 570), cv::FONT_HERSHEY_TRIPLEX, 0.5, cv::Scalar(0, 0, 0), 1, CV_AA);
	cv::putText(_img, "[Amplitude]", cv::Point(20 , 20), cv::FONT_HERSHEY_TRIPLEX, 0.5, cv::Scalar(0, 0, 0), 1, CV_AA);

	// x
	cv::putText(_img, "100", cv::Point(startX + (2-2) * timesX, 570), cv::FONT_HERSHEY_TRIPLEX, 0.5, cv::Scalar(0, 0, 0), 1, CV_AA);
	cv::putText(_img, "200", cv::Point(startX + (log10(200) - 2) * timesX, 570), cv::FONT_HERSHEY_TRIPLEX, 0.5, cv::Scalar(0, 0, 0), 1, CV_AA);
	cv::putText(_img, "500", cv::Point(startX + (log10(500)-2) * timesX, 570), cv::FONT_HERSHEY_TRIPLEX, 0.5, cv::Scalar(0, 0, 0), 1, CV_AA);
	cv::putText(_img, "800", cv::Point(startX + (log10(800) - 2) * timesX, 570), cv::FONT_HERSHEY_TRIPLEX, 0.5, cv::Scalar(0, 0, 0), 1, CV_AA);
	cv::putText(_img, "1000", cv::Point(startX + (3-2) * timesX, 570), cv::FONT_HERSHEY_TRIPLEX, 0.5, cv::Scalar(0, 0, 0), 1, CV_AA);
	cv::putText(_img, "10000", cv::Point(startX + (4 - 2) * timesX, 570), cv::FONT_HERSHEY_TRIPLEX, 0.5, cv::Scalar(0, 0, 0), 1, CV_AA);

	// y
	cv::putText(_img, "0.1", cv::Point(15, 550 - 1 * timesY), cv::FONT_HERSHEY_TRIPLEX, 0.5, cv::Scalar(0, 0, 0), 1, CV_AA);
	cv::putText(_img, "1", cv::Point(15, 550 - 2 * timesY), cv::FONT_HERSHEY_TRIPLEX, 0.5, cv::Scalar(0, 0, 0), 1, CV_AA);
	cv::putText(_img, "10", cv::Point(15, 550 - 3 * timesY), cv::FONT_HERSHEY_TRIPLEX, 0.5, cv::Scalar(0, 0, 0), 1, CV_AA);
}

void makeWav(const int xValue) 
{
	// Make sin wave
	waveSize = 1;
	double newFrequency = 0.0;
	double stepOfTime = 44100;
	uint time = 5;
	uint numberOfData = time * stepOfTime;
	vector<vector<double>> newData;
	newData.resize(numberOfData);

	/* 周波数の決定 */
	newFrequency = (xValue - 50) / timesX + 2; 	// double noramlizedDataX = startX + (dataX[i] - 2) * timesX;
	newFrequency = pow(10, newFrequency);
	cout << newFrequency << endl;

	for (int i = 0; i < numberOfData; ++i) {
		newData[i].resize(1);
		newData[i][0] = (double)waveSize * (double)sin(2.0 * M_PI*newFrequency*(double)(i / stepOfTime));
		//cout << newData[i][0] << endl;
	}

	sl::Wav newWav;
	newWav.bits = 16;
	newWav.channel = 1;
	newWav.length = numberOfData;
	newWav.data = newData;
	newWav.fs = stepOfTime;

	string name = "./output/"+to_string(newFrequency) + "Hz.wav";
	sl::wavwrite(name, newWav);
}

int main()
{
	/* wavファイルの準備 */
#if BOLERO
	// Read wav file
	wav = sl::wavread("./input/6616_London_g-1.wav");
	sl::showWavData(wav);
#endif
#if SIN_WAVE
	// Make 440Hz sin wave
	waveSize = 1;
	uint frequency = 400;
	double stepOfTime = 44100;
	uint time = 1;
	vector<double> data440;
	for (int i = 0; i < time * stepOfTime; ++i) {
		data440.push_back(100 * sin(2 * M_PI*frequency*(double)(i / stepOfTime)));
	}
#endif

#if CF	/* 相関関数 */
	calculateAcf(wav.data);
	calculateCcf(wav.data);
#endif

	/* 周波数解析 */
	cv::Mat frequency;
	uint indexOfStart = 0;
	vector<double> dataYPrevious, cluesPrevious;
	double wPrevious = 0.0;

	// FFTW
	fftwf_complex *in, *out;
	fftwf_plan plan;


	// 音符ごとの読み込み
	while (indexOfStart < wav.length) 
	{
		bool isSingle = true;
		uint lengthOfSingle = wav.fs / 10; // 初期値を適当に選ぶ

		while (isSingle) 
		{
			vector<double> dataYCurrent, cluesCurrent;
			double wCurrent = 0.0;

#if BOLERO // 基本角周波数
			wCurrent = (double)wav.fs / (double)lengthOfSingle;
#endif
#if SIN_WAVE
			double w = 1.0 / time;
#endif

			// メモリ確保
			in = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*lengthOfSingle);
			out = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*lengthOfSingle);

			// 初期化
#if BOLERO
			for (int i = indexOfStart; i < (indexOfStart + lengthOfSingle); ++i) 
			{
				in[i][0] = wav.data[i][0];
				in[i][1] = 0.0;
				out[i][0] = 0.0;
				out[i][1] = 0.0;
			}
#endif
#if SIN_WAVE
			for (int i = 0; i < data440.size(); ++i) {
				in[i][0] = data440[i];
				in[i][1] = 0.0;
				out[i][0] = 0.0;
				out[i][1] = 0.0;
			}
#endif

			// １次元のフーリエ変換を実行
			plan = fftwf_plan_dft_1d(lengthOfSingle, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
			fftwf_execute(plan);
			fftwf_destroy_plan(plan);

			//std::cout << indexOfStart << " 番目から " << lengthOfSingle << " 個分" << endl;
			cout << "[時間] "<< (double)(indexOfStart + lengthOfSingle) / (double)wav.fs << " s"<< endl;

			for (int i = 0; i < (lengthOfSingle / 2); i++)
			{
				double re = out[i][0];         // 複素数の実数部
				double im = out[i][1];         // 複素数の虚数部
				double mag = sqrt(re*re + im*im);      // 大きさを計算

				dataYCurrent.push_back(mag);
			}

			/* 平滑化 */
			const int range = 10;

			for (int i = 0; i < dataYCurrent.size(); ++i) {
				double sum = 0.0;
				int count = 0;

				for (int j = 0; j < 2*range; ++j) {
					if ((i - range + j) > 0 && (i - range + j) < dataYCurrent.size()) {
						count++;
						sum += dataYCurrent[i - range + j];
					}
				}
				if (count != 0) {
					dataYCurrent[i] = (double)(sum / (double)count);
				}
			}


#if 1 // 0.1秒ごとに表示
			for (int i = 0; i < dataYPrevious.size(); ++i)
			{
				dataYCurrent[i] = log10(dataYCurrent[i] * 100000 * SCALE);
			}

			/* 周波数応答のグラフ描画 */
			vector<double> dataXCurrent;
			for (int i = 0; i < dataYCurrent.size(); ++i) {
				dataXCurrent.push_back(log10(i*wCurrent));
			}
			DrawFrequency(frequency, dataXCurrent, dataYCurrent);
			cv::imshow("frequency", frequency);
			cv::waitKey(0);
#endif
			// 解放
			fftwf_free(in);
			fftwf_free(out);

			if (isSingle) {
				// 更新
				dataYPrevious = dataYCurrent;
				cluesPrevious = cluesCurrent;
				wPrevious = wCurrent;
				lengthOfSingle += wav.fs / 10; // 0.1秒スキップ
			}
			else {
				indexOfStart += lengthOfSingle;
				break;
			}
		}

		/* 縦軸のスケール調整 */
#if BOLERO
		for (int i = 0; i < dataYPrevious.size(); ++i) 
		{
			dataYPrevious[i] = log10(dataYPrevious[i] * 100000 * SCALE);
		}
#endif
#if SIN_WAVE
		double max = 0.0;
		for (int i = 0; i < dataY.size(); ++i) {
			if (maxY < dataY[i])	maxY = dataY[i];
		}
		double scale = 100 / maxY;
		cout << scale << endl;
		for (int i = 0; i < dataY.size(); ++i) {
			dataY[i] = log10(dataY[i] * scale);
		}
#endif

		/* 周波数応答のグラフ描画 */
		vector<double> dataXPrevious;
		for (int i = 0; i < dataYPrevious.size(); ++i) {
			dataXPrevious.push_back(log10(i*wPrevious));
		}
		DrawFrequency(frequency, dataXPrevious, dataYPrevious);
		cv::imshow("frequency", frequency);

		// 表示のためにストップ
		cv::waitKey(0);
		cout << "-----------------------------------------------------" << endl;
	}

	// コールバック関数
	mouseParam mouseEvent;
	cv::setMouseCallback("frequency", CallBackMouseFunc, &mouseEvent);

	while (1) {
		cv::waitKey(20);
		//左クリックがあったら表示
		if (mouseEvent.event == cv::EVENT_LBUTTONDOWN) {
			//クリック後のマウスの座標を出力
			std::cout << mouseEvent.x << " , " << mouseEvent.y << std::endl;
			makeWav(mouseEvent.x);
		}
		//右クリックがあったら終了
		else if (mouseEvent.event == cv::EVENT_RBUTTONDOWN) {
			break;
		}
	}

	std::cout << "終わり" << std::endl;
	std::system("pause");
}