#ifndef _INC_MCT_RAMPFILTER
#define _INC_MCT_RAMPFILTER


namespace mct
{

	enum FBPFilterWindow
	{
		FILTER_RAMLAK,			//< default FBP filter
		FILTER_SHEPPLOGAN,		//< Shepp-Logan
		FILTER_COSINE,			//< Cosine
		FILTER_HAMMING,			//< Hamming filter
		FILTER_HANN,			//< Hann filter
		FILTER_G_HAMMING,		//< Hamming filter
		FILTER_TUKEY,			//< Tukey filter
		FILTER_LANCZOS,			//< Lanczos filter
		FILTER_TRIANGULAR,		//< Triangular filter
		FILTER_GAUSSIAN,		//< Gaussian filter
		FILTER_BARTLETTHANN,	//< Bartlett-Hann filter
		FILTER_BLACKMAN,		//< Blackman filter
		FILTER_NUTTALL,			//< Nuttall filter, continuous first derivative
		FILTER_BLACKMANHARRIS,	//< Blackman-Harris filter
		FILTER_BLACKMANNUTTALL,	//< Blackman-Nuttall filter
		FILTER_FLATTOP,			//< Flat top filter
		FILTER_KAISER,			//< Kaiser filter
		FILTER_PARZEN,			//< Parzen filter
	};

	enum FilterGeneration
	{
		FILTER_GENERATION_DIRECT,
		FILTER_GENERATION_INVERSE_FOURIER,
	};

	enum SignalDomain
	{
		DOMAIN_FREQUENCY,
		DOMAIN_SPATIAL,
	};

	//关于重建滤波核
	//CT重建中可根据噪声、空间分辨率选择不同的滤波核。用户可根据需要选择标准、高分辨率、低噪声等滤波核。
	//滤波核对噪声、空间分辨率有影响，空间分辨率越高噪声越高，空间分辨率低噪声低。Tradeoff Spatial resolution <-> Noise
	//重建滤波核的构建通常根据窗函数、截止频率等进行设计
}

#endif