/// @ref core
/// @file glm/ext/vector_double3.hpp

#pragma once
#include "../detail/type_vec3.hpp"

namespace glm
{
	/// @addtogroup core_vector
	/// @{

	/// 3 components vector of double-precision floating-point numbers.
	///
	/// @see <a href="http://www.opengl.org/registry/doc/GLSLangSpec.4.20.8.pdf">GLSL 4.20.8 specification, section 4.1.5 Vectors</a>
	typedef vec<3, double, defaultp>		dvec3;

#if GLM_CONFIG_SIMD == GLM_ENABLE
	typedef vec<3, double, unaligned_simd>		usimd_dvec3;
#endif 

	/// @}
}//namespace glm
