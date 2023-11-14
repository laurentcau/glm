/// @ref core
/// @file glm/ext/vector_double4.hpp

#pragma once
#include "../detail/type_vec4.hpp"

namespace glm
{
	/// @addtogroup core_vector
	/// @{

	/// 4 components vector of double-precision floating-point numbers.
	///
	/// @see <a href="http://www.opengl.org/registry/doc/GLSLangSpec.4.20.8.pdf">GLSL 4.20.8 specification, section 4.1.5 Vectors</a>
	typedef vec<4, double, defaultp>		dvec4;

#if GLM_CONFIG_SIMD == GLM_ENABLE
	typedef vec<4, double, unaligned_simd>	usimd_dvec4;
#endif 
	/// @}
}//namespace glm
