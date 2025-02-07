/// @ref core
/// @file glm/ext/vector_float3.hpp

#pragma once
#include "../detail/type_vec3.hpp"

namespace glm
{
	/// @addtogroup core_vector
	/// @{

	/// 3 components vector of single-precision floating-point numbers.
	///
	/// @see <a href="http://www.opengl.org/registry/doc/GLSLangSpec.4.20.8.pdf">GLSL 4.20.8 specification, section 4.1.5 Vectors</a>
	typedef vec<3, float, defaultp>		vec3;

	typedef vec<3, float, packed_highp>		packed_vec3;

#if GLM_CONFIG_SIMD == GLM_ENABLE
	typedef vec<3, float, unaligned_simd>	usimd_vec3;
#endif 

	/// @}
}//namespace glm
