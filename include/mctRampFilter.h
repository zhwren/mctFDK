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

	//�����ؽ��˲���
	//CT�ؽ��пɸ����������ռ�ֱ���ѡ��ͬ���˲��ˡ��û��ɸ�����Ҫѡ���׼���߷ֱ��ʡ����������˲��ˡ�
	//�˲��˶��������ռ�ֱ�����Ӱ�죬�ռ�ֱ���Խ������Խ�ߣ��ռ�ֱ��ʵ������͡�Tradeoff Spatial resolution <-> Noise
	//�ؽ��˲��˵Ĺ���ͨ�����ݴ���������ֹƵ�ʵȽ������
}

#endif