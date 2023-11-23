#pragma once

#include "setup.hpp"

namespace glm
{
	/// Qualify GLM types in term of alignment (packed, aligned) and precision in term of ULPs (lowp, mediump, highp)
	enum qualifier
	{
		packed_highp, ///< Typed data is tightly packed in memory and operations are executed with high precision in term of ULPs
		packed_mediump, ///< Typed data is tightly packed in memory  and operations are executed with medium precision in term of ULPs for higher performance
		packed_lowp, ///< Typed data is tightly packed in memory  and operations are executed with low precision in term of ULPs to maximize performance

#		if GLM_CONFIG_SIMD == GLM_ENABLE
			unaligned_simd_highp, ///< Typed data is unaligned SIMD optimizations and operations are executed with high precision in term of ULPs
			unaligned_simd_mediump, ///< Typed data is unaligned SIMD optimizations and operations are executed with high precision in term of ULPs for higher performance
			unaligned_simd_lowp, // ///< Typed data is unaligned SIMD optimizations and operations are executed with high precision in term of ULPs to maximize performance
#		endif

#		if GLM_CONFIG_ALIGNED_GENTYPES == GLM_ENABLE
			aligned_highp, ///< Typed data is aligned in memory allowing SIMD optimizations and operations are executed with high precision in term of ULPs
			aligned_mediump, ///< Typed data is aligned in memory allowing SIMD optimizations and operations are executed with high precision in term of ULPs for higher performance
			aligned_lowp, // ///< Typed data is aligned in memory allowing SIMD optimizations and operations are executed with high precision in term of ULPs to maximize performance
#		endif

		highp = packed_highp, ///< By default highp qualifier is also packed
		mediump = packed_mediump, ///< By default mediump qualifier is also packed
		lowp = packed_lowp, ///< By default lowp qualifier is also packed
		packed = packed_highp, ///< By default packed qualifier is also high precision

#		if GLM_CONFIG_SIMD == GLM_ENABLE
			unaligned_simd = unaligned_simd_highp, ///< By default unaligned_simd qualifier is also high precision
#		endif

#		if GLM_CONFIG_ALIGNED_GENTYPES == GLM_ENABLE
			aligned = aligned_highp, ///< By default aligned qualifier is also high precision
#		endif

#		if GLM_CONFIG_ALIGNED_GENTYPES == GLM_ENABLE && defined(GLM_FORCE_DEFAULT_ALIGNED_GENTYPES)
			defaultp = aligned_highp
#		else
#			if GLM_CONFIG_SIMD == GLM_ENABLE
				defaultp = unaligned_simd_highp
#			else
				defaultp = highp
#			endif
#		endif
		
	};

	typedef qualifier precision;

	template<length_t L, typename T, qualifier Q = defaultp> struct vec;
	template<length_t C, length_t R, typename T, qualifier Q = defaultp> struct mat;
	template<typename T, qualifier Q = defaultp> struct qua;

#	if GLM_HAS_TEMPLATE_ALIASES
		template <typename T, qualifier Q = defaultp> using tvec1 = vec<1, T, Q>;
		template <typename T, qualifier Q = defaultp> using tvec2 = vec<2, T, Q>;
		template <typename T, qualifier Q = defaultp> using tvec3 = vec<3, T, Q>;
		template <typename T, qualifier Q = defaultp> using tvec4 = vec<4, T, Q>;
		template <typename T, qualifier Q = defaultp> using tmat2x2 = mat<2, 2, T, Q>;
		template <typename T, qualifier Q = defaultp> using tmat2x3 = mat<2, 3, T, Q>;
		template <typename T, qualifier Q = defaultp> using tmat2x4 = mat<2, 4, T, Q>;
		template <typename T, qualifier Q = defaultp> using tmat3x2 = mat<3, 2, T, Q>;
		template <typename T, qualifier Q = defaultp> using tmat3x3 = mat<3, 3, T, Q>;
		template <typename T, qualifier Q = defaultp> using tmat3x4 = mat<3, 4, T, Q>;
		template <typename T, qualifier Q = defaultp> using tmat4x2 = mat<4, 2, T, Q>;
		template <typename T, qualifier Q = defaultp> using tmat4x3 = mat<4, 3, T, Q>;
		template <typename T, qualifier Q = defaultp> using tmat4x4 = mat<4, 4, T, Q>;
		template <typename T, qualifier Q = defaultp> using tquat = qua<T, Q>;
#	endif

namespace detail
{
	template<glm::qualifier P>
	struct is_aligned
	{
		static const bool value = false;
	};

#	if GLM_CONFIG_ALIGNED_GENTYPES == GLM_ENABLE
		template<>
		struct is_aligned<glm::aligned_lowp>
		{
			static const bool value = true;
		};

		template<>
		struct is_aligned<glm::aligned_mediump>
		{
			static const bool value = true;
		};

		template<>
		struct is_aligned<glm::aligned_highp>
		{
			static const bool value = true;
		};
#	endif

		template<glm::qualifier P>
		struct use_simd
		{
			static const bool value = false;
		};

		// convert unaligned simd type to aligned simd type
		template<qualifier Q>
		struct to_aligned {};

#if GLM_CONFIG_SIMD == GLM_ENABLE
		template<>
		struct use_simd<glm::unaligned_simd_lowp>
		{
			static const bool value = true;
		};

		template<>
		struct use_simd<glm::unaligned_simd_mediump>
		{
			static const bool value = true;
		};

		template<>
		struct use_simd<glm::unaligned_simd_highp>
		{
			static const bool value = true;
		};

		template<>
		struct use_simd<glm::aligned_lowp>
		{
			static const bool value = true;
		};

		template<>
		struct use_simd<glm::aligned_mediump>
		{
			static const bool value = true;
		};

		template<>
		struct use_simd<glm::aligned_highp>
		{
			static const bool value = true;
		};

		template<>
		struct to_aligned<unaligned_simd_highp>
		{
			enum { value = aligned_highp };
		};

		template<>
		struct to_aligned<unaligned_simd_lowp>
		{
			enum { value = aligned_lowp };
		};

		template<>
		struct to_aligned<unaligned_simd_mediump>
		{
			enum { value = aligned_mediump };
		};

		template<>
		struct to_aligned<aligned_highp>
		{
			enum { value = aligned_highp };
		};

		template<>
		struct to_aligned<aligned_lowp>
		{
			enum { value = aligned_lowp };
		};

		template<>
		struct to_aligned<aligned_mediump>
		{
			enum { value = aligned_mediump };
		};

#endif

	template<length_t L, typename T, bool is_aligned, bool use_simd=false >
	struct storage
	{
		typedef struct type {
			T data[L];
		} type;
	};

#	if GLM_HAS_ALIGNOF
		template<length_t L, typename T>
		struct storage<L, T, true>
		{
			typedef struct alignas(L * sizeof(T)) type {
				T data[L];
			} type;
		};

		template<typename T>
		struct storage<3, T, true>
		{
			typedef struct alignas(4 * sizeof(T)) type {
				T data[4];
			} type;
		};

		template<>
		struct storage<3, bool, true, true>
		{
			typedef struct alignas(4 * sizeof(bool)) type {
				bool data[4];
			} type;
		};
#	endif

#	if GLM_ARCH & GLM_ARCH_SSE2_BIT
	template<>
	struct storage<4, float, true, true>
	{
		typedef glm_f32vec4 type;
	};

	template<>
	struct storage<4, float, false, true>
	{
		typedef struct type{
			float data[4];
			GLM_DEFAULTED_DEFAULT_CTOR_QUALIFIER GLM_CONSTEXPR type() GLM_DEFAULT;
			GLM_FUNC_QUALIFIER type(glm_f32vec4 v){_mm_storeu_ps(data, v);}
			GLM_FUNC_QUALIFIER operator glm_f32vec4() const {return _mm_loadu_ps(data);}
		} type;
	};


	template<>
	struct storage<4, int, true, true>
	{
		typedef glm_i32vec4 type;
	};

	template<>
	struct storage<4, int, false, true>
	{
		struct type
		{
			int data[4];
			GLM_DEFAULTED_DEFAULT_CTOR_QUALIFIER GLM_CONSTEXPR type() GLM_DEFAULT;
			GLM_FUNC_QUALIFIER type(glm_i32vec4 v) {
				_mm_storeu_si128((__m128i*)data, v); 
			}
			GLM_FUNC_QUALIFIER operator glm_i32vec4() const {
				return _mm_loadu_si128((__m128i*)data); 
			}
		};
	};

	template<>
	struct storage<4, unsigned int, true, true>
	{
		typedef glm_u32vec4 type;
	};

	template<>
	struct storage<4, unsigned int, false, true>
	{
		struct type
		{
			unsigned int data[4];
			GLM_DEFAULTED_DEFAULT_CTOR_QUALIFIER GLM_CONSTEXPR type() GLM_DEFAULT;
			GLM_FUNC_QUALIFIER type(glm_i32vec4 v) { _mm_storeu_si128((__m128i*)data, v); }
			GLM_FUNC_QUALIFIER operator glm_i32vec4() const { return _mm_loadu_si128((__m128i*)data); }
		};
	};

	template<>
	struct storage<3, float, true, true>
	{
		struct type
		{
			glm_f32vec4 data;
			//GLM_DEFAULTED_DEFAULT_CTOR_QUALIFIER GLM_CONSTEXPR type() { data = _mm_setzero_ps(); };
			GLM_DEFAULTED_DEFAULT_CTOR_QUALIFIER GLM_CONSTEXPR type() GLM_DEFAULT;
			GLM_FUNC_QUALIFIER type(glm_f32vec4 v) { data = v; }
			GLM_FUNC_QUALIFIER operator glm_f32vec4() const { return data; }
		};

	};

	template<>
	struct storage<3, float, false, true>
	{
		typedef struct type {
			float data[3];
			GLM_DEFAULTED_DEFAULT_CTOR_QUALIFIER GLM_CONSTEXPR type() GLM_DEFAULT;
			GLM_FUNC_QUALIFIER type(glm_f32vec4 v) {

#ifdef _MSC_VER 
				data[0] = v.m128_f32[0];
				data[1] = v.m128_f32[1];
				data[2] = v.m128_f32[2];
#else
				// other way to do it but generate more instructions
				_mm_store_sd(reinterpret_cast<double*>(data), _mm_castps_pd(v));
				__m128 z = _mm_shuffle_ps(v, v, _MM_SHUFFLE(2, 2, 2, 2));
				_mm_store_ss(&data[2], z);
#endif
				// other way to do it but generate more instructions
				//__m128 T1 = _mm_shuffle_ps(v, v, _MM_SHUFFLE(1, 1, 1, 1));
				//__m128 T2 = _mm_shuffle_ps(v, v, _MM_SHUFFLE(2, 2, 2, 2));
				//_mm_store_ss(&data[0], v);
				//_mm_store_ss(&data[1], T1);
				//_mm_store_ss(&data[2], T2);

			}
	
			GLM_FUNC_QUALIFIER operator glm_f32vec4() const {
				return _mm_set_ps(data[2], data[2], data[1], data[0]);

				// other way to do it but generate more instructions
				//__m128 xy = _mm_castpd_ps(_mm_load_sd(reinterpret_cast<const double*>(data)));
				//__m128 z = _mm_load_ss(&data[2]);
				//return _mm_movelh_ps(xy, z);
			}
		} type;
	};


	template<>
	struct storage<3, int, false, true>
	{
		typedef struct type {
			int data[3];
			GLM_DEFAULTED_DEFAULT_CTOR_QUALIFIER GLM_CONSTEXPR type() GLM_DEFAULT;
			GLM_FUNC_QUALIFIER type(glm_i32vec4 v) {
#ifdef _MSC_VER 
				data[0] = v.m128i_i32[0];
				data[1] = v.m128i_i32[1];
				data[2] = v.m128i_i32[2];
#else
				_mm_store_sd(reinterpret_cast<double*>(data), _mm_castsi128_pd(v));
				__m128 z = _mm_shuffle_ps(_mm_castsi128_ps(v), _mm_castsi128_ps(v), _MM_SHUFFLE(2, 2, 2, 2));
				_mm_store_ss(reinterpret_cast<float*>(&data[2]), z);
#endif
			}
			GLM_FUNC_QUALIFIER operator glm_i32vec4() const {
#ifdef _MSC_VER 
				glm_i32vec4 v;
				v.m128i_i32[0] = data[0];
				v.m128i_i32[1] = data[1];
				v.m128i_i32[2] = data[2];
				v.m128i_i32[3] = 0;
				return v;
#else
				__m128 x = _mm_load_ss(reinterpret_cast<const float*>(&data[0]));
				__m128 y = _mm_load_ss(reinterpret_cast<const float*>(&data[1]));
				__m128 z = _mm_load_ss(reinterpret_cast<const float*>(&data[2]));
				__m128 xy = _mm_unpacklo_ps(x, y);
				return _mm_castps_si128(_mm_movelh_ps(xy, z));
#endif
			}
		} type;
	};

	template<>
	struct storage<3, unsigned int, false, true>
	{
		typedef struct type {
			unsigned int data[3];
			GLM_DEFAULTED_DEFAULT_CTOR_QUALIFIER GLM_CONSTEXPR type() GLM_DEFAULT;
			GLM_FUNC_QUALIFIER type(glm_i32vec4 v) {
#ifdef _MSC_VER 
				data[0] = v.m128i_i32[0];
				data[1] = v.m128i_i32[1];
				data[2] = v.m128i_i32[2];
#else
				_mm_store_sd(reinterpret_cast<double*>(data), _mm_castsi128_pd(v));
				__m128 z = _mm_shuffle_ps(_mm_castsi128_ps(v), _mm_castsi128_ps(v), _MM_SHUFFLE(2, 2, 2, 2));
				_mm_store_ss(reinterpret_cast<float*>(&data[2]), z);
#endif
			}
			GLM_FUNC_QUALIFIER operator glm_i32vec4() const {
#ifdef _MSC_VER 
				glm_i32vec4 v;
				v.m128i_i32[0] = data[0];
				v.m128i_i32[1] = data[1];
				v.m128i_i32[2] = data[2];
				v.m128i_i32[3] = 0;
				return v;
#else
				__m128 x = _mm_load_ss(reinterpret_cast<const float*>(&data[0]));
				__m128 y = _mm_load_ss(reinterpret_cast<const float*>(&data[1]));
				__m128 z = _mm_load_ss(reinterpret_cast<const float*>(&data[2]));
				__m128 xy = _mm_unpacklo_ps(x, y);
				return _mm_castps_si128(_mm_movelh_ps(xy, z));
#endif
			}
		} type;
	};

	template<>
	struct storage<3, int, true, true>
	{
		typedef glm_i32vec4 type;
	};

	template<>
	struct storage<3, unsigned int, true, true>
	{
		typedef glm_i32vec4 type;
	};

	template<>
	struct storage<2, double, true, true>
	{
		typedef glm_f64vec2 type;
	};

	template<>
	struct storage<2, double, false, true>
	{
		struct type
		{
			double data[2];
			GLM_DEFAULTED_DEFAULT_CTOR_QUALIFIER GLM_CONSTEXPR type() GLM_DEFAULT;
			type(glm_f64vec2 v) { _mm_storeu_pd(data, v); }
			operator glm_f64vec2() const { return _mm_loadu_pd(data); }
		};
	};

	template<>
	struct storage<2, detail::int64, true, true>
	{
		typedef glm_i64vec2 type;
	};

	template<>
	struct storage<2, detail::uint64, true, true>
	{
		typedef glm_u64vec2 type;
	};


	template<>
	struct storage<3, detail::uint64, true, true>
	{
		typedef glm_u64vec2 type;
	};

	template<>
	struct storage<4, double, false, true>
	{
		struct type
		{
			double data[4];
			GLM_DEFAULTED_DEFAULT_CTOR_QUALIFIER GLM_CONSTEXPR type() GLM_DEFAULT;
#	if (GLM_ARCH & GLM_ARCH_AVX_BIT)
			type(glm_f64vec4 v) { _mm256_storeu_pd(data, v); }
			operator glm_f64vec4() const { return _mm256_loadu_pd(data); }
#else
			GLM_FUNC_QUALIFIER GLM_CONSTEXPR glm_f64vec2 getv(int i) const {
				return _mm_loadu_pd(i == 0 ? data : &data[2]);
			}
			GLM_FUNC_QUALIFIER GLM_CONSTEXPR void setv(int i, const glm_f64vec2& v) {
				_mm_storeu_pd(i == 0 ? data : &data[2], v);
			}
#endif
		};
	};


	template<>
	struct storage<4, double, true, true>
	{
#	if (GLM_ARCH & GLM_ARCH_AVX_BIT)
		typedef glm_f64vec4 type;
#	else
		struct type
		{
			glm_f64vec2 data[2];
			GLM_CONSTEXPR glm_f64vec2 getv(int i) const {
				return data[i];
			}
			GLM_CONSTEXPR void setv(int i, const glm_f64vec2& v) {
				data[i] = v;
			}
		};
#	endif
	};


	template<>
	struct storage<3, double, true, true> : public storage<4, double, true, true>
	{};

	template<>
	struct storage<3, double, false, true>
	{
		struct type
		{
			double data[3];
			GLM_DEFAULTED_DEFAULT_CTOR_QUALIFIER GLM_CONSTEXPR type() GLM_DEFAULT;
#	if (GLM_ARCH & GLM_ARCH_AVX_BIT)
			GLM_FUNC_QUALIFIER type(glm_f64vec4 v) {
				__m256d T1 = _mm256_permute_pd(v, 1);
				_mm_store_sd(&data[0], _mm256_castpd256_pd128(v));
				_mm_store_sd(&data[1], _mm256_castpd256_pd128(T1));
				_mm_store_sd(&data[2], _mm256_extractf128_pd(v, 1));
			}
			GLM_FUNC_QUALIFIER operator glm_f64vec4() const {
				return _mm256_set_pd(data[2], data[2], data[1], data[0]);
			}
#	else
			GLM_CONSTEXPR glm_f64vec2 getv(int i) const {
				if (i == 0) 
					return _mm_loadu_pd(data); 
				else
					return _mm_load_sd(&data[2]);
			}
			GLM_CONSTEXPR void setv(int i, const glm_f64vec2 &v)
			{
				if (i == 0)
					_mm_storeu_pd(data, v);
				else
					_mm_store_sd(&data[2], v);
			}
#endif
		};

	};

#	endif


#	if (GLM_ARCH & GLM_ARCH_AVX2_BIT)
	template<>
	struct storage<4, detail::int64, true, true>
	{
		typedef glm_i64vec4 type;
	};

	template<>
	struct storage<4, detail::uint64, true, true>
	{
		typedef glm_u64vec4 type;
	};
#	endif

#	if GLM_ARCH & GLM_ARCH_NEON_BIT
	template<>
	struct storage<4, float, true, true>
	{
		typedef glm_f32vec4 type;
	};

	template<>
	struct storage<4, float, false, true>
	{
		typedef struct type {
			float data[4];
			GLM_DEFAULTED_DEFAULT_CTOR_QUALIFIER GLM_CONSTEXPR type() GLM_DEFAULT;
			inline type(glm_f32vec4 v) { vst1q_f32(reinterpret_cast<float*>(data), v); }
			inline operator glm_f32vec4() const { return vld1q_f32(reinterpret_cast<const float*>(data)); }
		} type;
	};

	template<>
	struct storage<3, float, true, true> : public storage<4, float, true, true>
	{};

	template<>
	struct storage<3, float, false, true>
	{
		typedef struct type {
			float data[3];
			GLM_DEFAULTED_DEFAULT_CTOR_QUALIFIER GLM_CONSTEXPR type() GLM_DEFAULT;
			inline type(glm_f32vec4 v) { 
				data[0] = vgetq_lane_f32(v, 0);
				data[1] = vgetq_lane_f32(v, 1);
				data[2] = vgetq_lane_f32(v, 2);
			}
			inline operator glm_f32vec4() const { 
				const float init[4] = {data[2],data[2],data[1],data[0]};
				return vld1q_f32(init);
			}
		} type;

	};

	template<>
	struct storage<4, int, true, true>
	{
		typedef glm_i32vec4 type;
	};

	template<>
	struct storage<4, int, false, true>
	{
		struct type
		{
			int data[4];
			GLM_DEFAULTED_DEFAULT_CTOR_QUALIFIER GLM_CONSTEXPR type() GLM_DEFAULT;
			type(glm_i32vec4 v) { vst1q_s32(data, v); }
			operator glm_i32vec4() const { return vld1q_s32(data); }
		};
	};

	template<>
	struct storage<3, int, true, true> : public storage<4, int, true, true>
	{};

	template<>
	struct storage<3, int, false, true>
	{
		struct type
		{
			int data[3];
			GLM_DEFAULTED_DEFAULT_CTOR_QUALIFIER GLM_CONSTEXPR type() GLM_DEFAULT;
			inline type(glm_i32vec4 v) {
				data[0] = vgetq_lane_s32(v, 0);
				data[1] = vgetq_lane_s32(v, 1);
				data[2] = vgetq_lane_s32(v, 2);
			}
			inline operator glm_i32vec4() const {
				const int init[4] = { data[2], data[2], data[1], data[0] };
				return vld1q_s32(init);
			}
		};
	};

	template<>
	struct storage<4, unsigned int, true, true>
	{
		typedef glm_u32vec4 type;
	};

	template<>
	struct storage<4, unsigned int, false, true>
	{
		struct type
		{
			unsigned int data[4];
			GLM_DEFAULTED_DEFAULT_CTOR_QUALIFIER GLM_CONSTEXPR type() GLM_DEFAULT;
			type(glm_u32vec4 v) { vst1q_u32(data, v); }
			operator glm_u32vec4() const { return vld1q_u32(data); }
		};
	};

	template<>
	struct storage<3, unsigned int, true, true> : public storage<4, unsigned int, true, true>
	{};

	template<>
	struct storage<3, unsigned int, false, true>
	{
		struct type
		{
			unsigned int data[3];
			GLM_DEFAULTED_DEFAULT_CTOR_QUALIFIER GLM_CONSTEXPR type() GLM_DEFAULT;
			inline type(glm_u32vec4 v) {
				data[0] = vgetq_lane_u32(v, 0);
				data[1] = vgetq_lane_u32(v, 1);
				data[2] = vgetq_lane_u32(v, 2);
			}
			inline operator glm_u32vec4() const {
				const unsigned int init[4] = { data[2], data[2], data[1], data[0] };
				return vld1q_u32(init);
			}
		};
	};

	template<>
	struct storage<3, double, true, true>
	{
		typedef struct alignas(4 * sizeof(double)) type {
			double data[4];
		} type;
	};

#	endif

	enum genTypeEnum
	{
		GENTYPE_VEC,
		GENTYPE_MAT,
		GENTYPE_QUAT
	};

	template <typename genType>
	struct genTypeTrait
	{};

	template <length_t C, length_t R, typename T>
	struct genTypeTrait<mat<C, R, T> >
	{
		static const genTypeEnum GENTYPE = GENTYPE_MAT;
	};

	template<typename genType, genTypeEnum type>
	struct init_gentype
	{
	};

	template<typename genType>
	struct init_gentype<genType, GENTYPE_QUAT>
	{
		GLM_FUNC_QUALIFIER GLM_CONSTEXPR static genType identity()
		{
			return genType(1, 0, 0, 0);
		}
	};

	template<typename genType>
	struct init_gentype<genType, GENTYPE_MAT>
	{
		GLM_FUNC_QUALIFIER GLM_CONSTEXPR static genType identity()
		{
			return genType(1);
		}
	};
}//namespace detail
}//namespace glm
