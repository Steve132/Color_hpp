#ifndef COLORSPACES_HPP
#define COLORSPACES_HPP

#include<array>
#include<tuple>
#include<utility>
#include<type_traits>
#include<numeric>
#include<ratio>
#include<cmath>

namespace {
template<class T>
static inline T sconv(T x){ return x; }

/*template<class T,class U,
		 std::enable_if_t<std::is_integral_v<T> && std::is_integral_v<U>>
		 >
static inline T sconv(U x) {
	static constexpr std::ratio<std::numeric_limits<T>::max(),std::ratio<std::numeric_limits<U>::max()> scalar;

	return x
}*/
static inline float sconv(uint8_t x){ return static_cast<float>(x)/255.0f; }
static inline uint8_t sconv(float x){ return static_cast<uint8_t>(x*255.0f); }
}

template<class Dst,class Src,class ...Args>
struct color_convert
{};

template<class T>
struct sRGB{
	T R,G,B;

	sRGB()=default;
	sRGB(T tR,T tG,T tB):R(tR),G(tG),B(tB){}

	template<class U>
	explicit sRGB(const sRGB<U>& o):sRGB{sconv(o.R),sconv(o.G),sconv(o.B)}
	{}
	operator std::array<T,3>() const
	{return std::array<T,3>{R,G,B};}

	template<class Other,class ...Args>
	explicit sRGB(const Other& o,Args&&...args):sRGB(color_convert<sRGB,Other,Args...>::convert(o,std::forward(args)...))
	{}
};
inline sRGB<float> cube_clamp(sRGB<float> f)
{
	float mx=std::max(std::max(f.R,f.G),f.B);
	if(mx > 1.0f)
	{
		mx=1.0f/mx;
		f.R*=mx;f.G*=mx;f.B*=mx;
	}
	return sRGB<float>{
		std::max(0.0f,f.R),
		std::max(0.0f,f.G),
		std::max(0.0f,f.B)
	};
}

template<class T>
struct CIEXYZ{
	T X,Y,Z;
	CIEXYZ()=default;
	CIEXYZ(T tX,T tY,T tZ):X(tX),Y(tY),Z(tZ){}
	template<class U>
	explicit CIEXYZ(const CIEXYZ<U>& o):CIEXYZ{sconv(o.X),sconv(o.Y),sconv(o.Z)}
	{}
	operator std::array<T,3>() const
	{return std::array<T,3>{X,Y,Z};}

	template<class Other,class ...Args>
	explicit CIEXYZ(const Other& o,Args&&...args):CIEXYZ(color_convert<CIEXYZ,Other,Args...>::convert(o,std::forward(args)...))
	{}
};


template<class T>
struct CIELab{
	T L,a,b;

	CIELab()=default;
	CIELab(T tL,T ta,T tb):L(tL),a(ta),b(tb){}

	template<class U>
	explicit CIELab(const CIELab<U>& o):CIELab{sconv(o.R),sconv(o.G),sconv(o.B)}
	{}
	operator std::array<T,3>() const
	{return std::array<T,3>{L,a,b};}

	template<class Other,class ...Args>
	explicit CIELab(const Other& o,Args&&...args):CIELab(color_convert<CIELab,Other,Args...>::convert(o,std::forward(args)...))
	{}
};

template<size_t N,class Func>
static constexpr auto lutgen(Func&& f)->
	std::array<decltype(f(size_t{})),N>
{
	std::array<decltype(f(size_t{})),N> lut;
	for(size_t i=0;i<N;i++) lut[i]=f(i);
	return lut;
}

namespace{

template<unsigned int N>
static inline float polyval(float x,const std::array<float,N>& coeffs){
	float cx=1.0f;
	return std::accumulate(coeffs.begin()+1,coeffs.end(),
	coeffs[0],[&cx,x](float v,float c)
	{
		cx*=x;
		return v+c*cx;
	}
	);
}
}

static constexpr float sRGBgamma(float x)
{
	if(x <= 0.04045f) return x/(12.92f);

	return ::polyval<4>(x,std::array<float,4>{     3.9345289214513454e-003f,
									-1.3128636574970765e-002f,
									 7.2637968830483990e-001f,
									 2.8351493886165202e-001f
					});
}
static constexpr float sRGBigamma(float x)
{
	if(x <= 0.0031308) return 12.92f*x;
	return ::polyval<6>(x,std::array<float,6>{
		1.2879089912942732e-001,
		2.6989743978570759e+000,
	   -6.2616991567126581e+000,
		1.0029861160531468e+001,
	   -8.2519091430162810e+000,
		2.6583005818050025e+000
	});
}

static constexpr float sRGBgamma(uint8_t x)
{
	std::array<float,256> lut=lutgen<256>([](size_t i){
			float iv=static_cast<float>(i)/255.0f;
			return sRGBgamma(iv);});
	return lut[x];
}



template<class T,class U>
struct color_convert<CIEXYZ<T>,sRGB<U>>
{
	static inline CIEXYZ<T> convert(const sRGB<U>& o)
	{
		float Rp=sRGBgamma(o.R);
		float Gp=sRGBgamma(o.G);
		float Bp=sRGBgamma(o.B);
		return CIEXYZ<float>{
						0.4124564f*Rp+0.3575761f*Gp+0.1804375f*Bp,
						0.2126729f*Rp+0.7151522f*Gp+0.0721750f*Bp,
						0.0193339f*Rp+0.1191920f*Gp+0.9503041f*Bp
		};
	}
};

template<class T,class U>
struct color_convert<sRGB<T>,CIEXYZ<U>>
{
	static inline sRGB<T> convert(const CIEXYZ<U>& o)
	{
		CIEXYZ<float> f=o;

		float Rp= 3.2405f*f.X-1.5372f*f.Y-0.4986f*f.Z;
		float Gp=-0.9689f*f.X+1.8758f*f.Y+0.0415f*f.Z;
		float Bp= 0.0557f*f.X-0.2040f*f.Y+1.0570f*f.Z;

		return sRGB<float>{sRGBigamma(Rp),sRGBigamma(Gp),sRGBigamma(Bp)};
	}
};





template<class T>
struct Lightness{
	T L; //Lum may be lightness, saturation, etc.

	Lightness()=default;
	Lightness(T tL):
		L(tL){}

	template<class U>
	explicit Lightness(const Lightness<U>& o):Lightness{sconv(o.L)}
	{}

	template<class Other,class ...Args>
	explicit Lightness(const Other& o,Args&&...args):
		Lightness(color_convert<Lightness,Other,Args...>::convert(o,std::forward(args)...))
	{}
};

template<class T>
struct Intensity{
	T I; //Lum may be lightness, saturation, etc.

	Intensity()=default;
	Intensity(T tI):
		I(tI){}

	template<class U>
	explicit Intensity(const Intensity<U>& o):Intensity{sconv(o.L)}
	{}

	template<class Other,class ...Args>
	explicit Intensity(const Other& o,Args&&...args):
		Intensity(color_convert<Intensity,Other,Args...>::convert(o,std::forward(args)...))
	{}
};

template<class T>
struct Luma{
	T Y; //Lum may be lightness, saturation, etc.

	Luma()=default;
	Luma(T tY):
		Y(tY){}

	template<class U>
	explicit Luma(const Luma<U>& o):Luma{sconv(o.Y)}
	{}

	template<class Other,class ...Args>
	explicit Luma(const Other& o,Args&&...args):
		Luma(color_convert<Luma,Other,Args...>::convert(o,std::forward(args)...))
	{}
};



template<class T,class U>
struct color_convert<Lightness<T>,sRGB<U>>
{
	static inline Lightness<T> convert(const sRGB<U>& o)
	{
		sRGB<float> f=o;
		float M=std::max(std::max(f.R,f.G),f.B);
		float m=std::min(std::min(f.R,f.G),f.B);
		return Lightness<float>{.5f*(M+m)};
	}
};


template<class T,class U>
struct color_convert<Intensity<T>,sRGB<U>>
{
	static inline Intensity<T> convert(const sRGB<U>& o)
	{
		sRGB<float> f=o;

		return Intensity<float>{(f.R+f.G+f.B)/3.0f};
	}
};


template<class T>
struct HueChroma{
	T H,C;

	HueChroma()=default;
	HueChroma(T tH,T tC):
		H(tH),C(tC){}

	template<class U>
	explicit HueChroma(const HueChroma<U>& o):HueChroma{sconv(o.H),sconv(o.C)}
	{}
	operator std::array<T,2>() const
	{return std::array<T,2>{H,C};}

	template<class Other,class ...Args>
	explicit HueChroma(const Other& o,Args&&...args):
		HueChroma(color_convert<HueChroma,Other,Args...>::convert(o,std::forward(args)...))
	{}
};



template<class T,class U>
struct color_convert<HueChroma<T>,sRGB<U>>
{
	static inline HueChroma<T> convert(const sRGB<U>& o)
	{
		sRGB<float> f=o;
		constexpr float PI=(float)M_PI;


		unsigned int selector= (static_cast<unsigned int>(f.R > f.G) << 2)
							| (static_cast<unsigned int>(f.G > f.B)  << 1)
							| (static_cast<unsigned int>(f.B > f.R)  << 0);
		float offset=0.0f;
		float H,C;

		switch(selector){
			case 0x0:
			case 0x7:	//0x7 is impossible
			{
				return HueChroma<float>{0.0f,0.0f};
			}
			case 0x4:
				offset=6.0f; //fallthrough with the 6 if g < b
			case 0x6:
			{
				C=f.R-f.B;    //(g < b ? 6 : 0)
				H=((f.G - f.B) / C) + offset; // justification: https://stackoverflow.com/a/39147465
				break;
			}
			case 0x2:
			case 0x3:
			{
				offset=2.0f;
				C=f.G-f.R;
				H=((f.B - f.R) / C) + offset;
				break;
			}
			case 0x1:
			case 0x5:
			{
				offset=4.0f;
				C=f.B-f.G;
				H=((f.R - f.G) / C) + offset;
			}
		};

		return HueChroma<float>{H*PI/3.0f,C};
	}
};

template<class T,class U>
struct color_convert<sRGB<T>,HueChroma<U>>
{
	static inline sRGB<T> convert(const HueChroma<U>& o,const Lightness<U>& l);
};


template<class T>
struct CircularHueChroma{
	T H,C; //Lum may be lightness, saturation, etc.

	CircularHueChroma()=default;
	CircularHueChroma(T tH,T tC):
		H(tH),C(tC){}

	template<class U>
	explicit CircularHueChroma(const CircularHueChroma<U>& o):CircularHueChroma{sconv(o.H),sconv(o.C)}
	{}
	operator std::array<T,2>() const
	{return std::array<T,2>{H,C};}

	template<class Other,class ...Args>
	explicit CircularHueChroma(const Other& o,Args&&...args):
		CircularHueChroma(color_convert<CircularHueChroma,Other,Args...>::convert(o,std::forward(args)...))
	{}


};

template<class T,class U>
struct color_convert<CircularHueChroma<T>,sRGB<U>>
{
	static inline CircularHueChroma<T> convert(const sRGB<U>& o)
	{
		sRGB<float> f=o;
		static constexpr float sin60=0.86602540378f;
		float alpha=f.R-0.5f*f.G-0.5f*f.B;
		float beta=(f.G-f.B)*sin60;
		float hue=std::atan2(beta,alpha);
		float chroma=std::sqrt(alpha*alpha+beta*beta);
		return CircularHueChroma<float>{hue,chroma};
	}
};
template<class T,class U>
struct color_convert<sRGB<T>,CircularHueChroma<U>>
{
	static inline sRGB<T> convert(const CircularHueChroma<U>& o,const Intensity<U>& l)
	{
		float H=sconv<float>(o.H);
		float C=sconv<float>(o.C);
		float I=sconv<float>(l.I);
		static constexpr float sin60=0.86602540378f;
		float alpha=C*cosf(H);
		float beta=C*sinf(H);
		float GmB=beta/sin60;
		float I3=3.0f*I;
		float R=(I3+2.0f*alpha)/2.0f;
		float GpB=-2.0f*(alpha-R);
		float G=GmB+GpB;
		float B=GpB-G;
		return sRGB<float>{R,G,B};
	}
};



static constexpr float lab_f(float x)
{
	if(x < 0.00885645167f) return x*7.78703703704f+0.13793103448f;
	return polyval<7>(x,std::array<float,7>{
		2.1455433364075635e-001f,
		3.5416223098428588e+000f,
	   -1.3602287189198677e+001f,
		3.4538715501222448e+001f,
	   -4.8821033387852289e+001f,
		3.5264766269008852e+001f,
	   -1.0138949423356147e+001f
	});
}
static constexpr float lab_if(float x)
{
	if(x < 0.20689655172f) return x*0.12841854934f-0.01771290335f;
	return x*x*x;
}

template<class T,class U>
struct color_convert<CIELab<T>,CIEXYZ<U>>
{
	static inline CIELab<T> convert(const CIEXYZ<U>& o,
			const CIEXYZ<float>& illum={0.950489f,1.0f,1.08884f} //Assume D65
		)
	{
		CIEXYZ<float> f{lab_f(o.X/illum.X),
						lab_f(o.Y/illum.Y),
					lab_f(o.Z/illum.Z)};

		float Lp=1.16f*f.X-0.16;
		float ap=5.0f*(f.X-f.Y);
		float bp=2.0f*(f.Y-f.Z);

		return CIELab<float>{Lp,ap,bp};
	}
};
template<class T,class U>
struct color_convert<CIEXYZ<U>,CIELab<T>>
{
	static inline CIEXYZ<U> convert(const CIELab<U>& o,
									const CIEXYZ<float>& illum={0.950489f,1.0f,1.08884f} //Assume D65
		)
	{
		CIELab<float> f{
			(o.L+0.16f)/1.16f,
			o.a/5.0f,
			o.b/2.0f
		};

		return CIEXYZ<float>{
			illum.X * lab_if(f.L + f.a),
			illum.Y * lab_if(f.L),
			illum.Z * lab_if(f.L - f.b)
		};
	}
};


//Todo: maybe all the circular models should be based on a single common base class?
template<class T>
struct CIELCh_ab{
	T L,C,h;

	CIELCh_ab()=default;
	CIELCh_ab(T tL,T tC,T th):
		L(tL),C(tC),h(th){}

	template<class U>
	explicit CIELCh_ab(const CIELCh_ab<U>& o):CIELCh_ab{sconv(o.L),sconv(o.C),sconv(o.h)}
	{}
	operator std::array<T,3>() const
	{return std::array<T,3>{L,C,h};}

	template<class Other,class ...Args>
	explicit CIELCh_ab(const Other& o,Args&&...args):
		CIELCh_ab(color_convert<CIELCh_ab,Other,Args...>::convert(o,std::forward(args)...))
	{}
};
template<class T,class U>
struct color_convert<CIELCh_ab<T>,CIELab<U>>
{
	static inline CIELCh_ab<T> convert(const CIELab<U>& o)
	{
		return CIELCh_ab<T>{o.L,std::hypot(o.a,o.b),std::atan2(o.b,o.a)+(float)M_PI};
	}
};
template<class T,class U>
struct color_convert<CIELab<T>,CIELCh_ab<U>>
{
	static inline CIELab<T> convert(const CIELCh_ab<U>& o)
	{
		return CIELab<T>{o.L,-o.C*std::cos(o.h),-o.C*std::sin(o.h)};
	}
};


using Color3f=std::array<float,3>;
using Color4f=std::array<float,4>;

#endif

/*
inline Color3f pix2rgT(const std::array<uint8_t,3>& RGB)
{
	Color3f rgb(RGB);
	float r=rgb[0],g=rgb[1],b=rgb[2];
	float T=r+g+b;
	if(T==0.0f) return Color3f{0.0f,0.0f,0.0f};
    else return Color3f{r/T,g/T,T/3.0f};
}
inline std::array<uint8_t,3> rgT2pix(const Color3f& rgT)
{
	float r=rgT[0],g=rgT[1];
    float T=3.0f*rgT[2];
    float b=T-r-g;
	Color3f rgb{r*T,g*T,b*T};
	return std::array<uint8_t,3>(rgb);
}


std::array<float,3> rgb2wbf(const std::array<uint8_t,3>& in)
{
	sRGB<float> sRGBin{sRGB<uint8_t>{in[0],in[1],in[2]}};
	CircularHueChroma<float> hc{sRGBin};
	Intensity<float> l{sRGBin};

	return std::array<float,3>{hc.H,hc.C,l.I};
}



template<typename T1, typename T2, typename F>
void cimg_color_transform(const CImg<T1>& input, CImg<T2>& output, F f)
{
	cimg_forXY(output,x,y)
	{
		std::array<T1,3> in_pixel{input(x,y,0),input(x,y,1),input(x,y,2)};
		std::array<T2,3> out_pixel = f(in_pixel);
		output(x,y,0) = out_pixel[0];
		output(x,y,1) = out_pixel[1];
		output(x,y,2) = out_pixel[2];
	}

}
*/
