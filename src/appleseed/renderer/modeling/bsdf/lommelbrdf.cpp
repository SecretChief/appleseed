
//
// This source file is part of appleseed.
// Visit http://appleseedhq.net/ for additional information and resources.
//
// This software is released under the MIT license.
//
// Copyright (c) 2014 Hans Hoogenboom, The appleseedhq Organization
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.
//

// Interface header.
#include "lommelbrdf.h"

// appleseed.renderer headers.
#include "renderer/modeling/bsdf/bsdf.h"
#include "renderer/modeling/bsdf/bsdfwrapper.h"

// appleseed.foundation headers.
#include "foundation/math/basis.h"
#include "foundation/math/sampling.h"
#include "foundation/math/vector.h"
#include "foundation/utility/containers/dictionary.h"
#include "foundation/utility/containers/specializedarrays.h"

// Standard headers.
#include <cmath>

// Forward declarations.
namespace foundation    { class AbortSwitch; }
namespace renderer      { class Assembly; }
namespace renderer      { class Project; }

using namespace foundation;
using namespace std;

namespace renderer
{

namespace
{
    //
    // Lommel BRDF.
    //

    const char* Model = "lommel_brdf";

    class LommelBRDFImpl
      : public BSDF
    {
      public:
        LommelBRDFImpl(
            const char*         name,
            const ParamArray&   params)
          : BSDF(name, Reflective, Diffuse, params)
        {
            m_inputs.declare("reflectance", InputFormatSpectralReflectance);
            m_inputs.declare("reflectance_multiplier", InputFormatScalar, "1.0");
            m_inputs.declare("thickness", InputFormatScalar, "0.1");
        }

        virtual void release() OVERRIDE
        {
            delete this;
        }

        virtual const char* get_model() const OVERRIDE
        {
            return Model;
        }

        void lommel_evaluate( 
            const double    c_in, 
            const double    c_on, 
            const double    reflectance_multiplier, 
            const double    thickness, 
            const Spectrum& reflectance,
            Spectrum&       value) const
        {
            const double a = 1.0 / (c_in + c_on);
            const double b = 1.0 - exp(-thickness * (1.0 / c_in + 1.0 / c_on));

            value = reflectance;
            value *= static_cast<float>(reflectance_multiplier * a * b * c_in);
        }

        FORCE_INLINE virtual Mode sample(
            SamplingContext&    sampling_context,
            const void*         data,
            const bool          adjoint,
            const bool          cosine_mult,
            const Vector3d&     geometric_normal,
            const Basis3d&      shading_basis,
            const Vector3d&     outgoing,
            Vector3d&           incoming,
            Spectrum&           value,
            double&             probability) const
        {
            // Compute the incoming direction in local space.
            sampling_context.split_in_place(2, 1);
            const Vector2d s = sampling_context.next_vector2<2>();
            const Vector3d wi = sample_hemisphere_cosine(s);

            // Transform the incoming direction to parent space.
            incoming = shading_basis.transform_to_parent(wi);

            // Compute the BRDF value.
            const InputValues* values = static_cast<const InputValues*>(data);
            const Vector3d& n   = shading_basis.get_normal();
            
            const double cos_in = dot(incoming, n);
            const double cos_on = dot(outgoing, n);
            
            if (cos_in == -cos_on || cos_in == 0.0 || cos_on <= 0.0)
                return Absorption;

            lommel_evaluate(cos_in, cos_on, values->m_thickness, values->m_reflectance_multiplier, values->m_reflectance, value);

            // Compute the probability density of the sampled direction.
            probability = wi.y * RcpPi;
            assert(probability > 0.0);

            // Return the scattering mode.
            return Diffuse;
        }

        FORCE_INLINE virtual double evaluate(
            const void*         data,
            const bool          adjoint,
            const bool          cosine_mult,
            const Vector3d&     geometric_normal,
            const Basis3d&      shading_basis,
            const Vector3d&     outgoing,
            const Vector3d&     incoming,
            const int           modes,
            Spectrum&           value) const
        {
            if (!(modes & Diffuse))
                return 0.0;

            // No reflection below the shading surface.
            const Vector3d& n = shading_basis.get_normal();
            const double cos_in = dot(incoming, n);
            if (cos_in < 0.0)
                return 0.0;

            // Compute the BRDF value.
            const InputValues* values = static_cast<const InputValues*>(data);
            const double cos_on = dot(outgoing, n);

            if (cos_in == -cos_on || cos_in == 0.0 || cos_on <= 0.0)
                return 0.0;

            lommel_evaluate(cos_in, cos_on, values->m_thickness, values->m_reflectance_multiplier, values->m_reflectance, value);

            // Return the probability density of the sampled direction.
            return cos_in * RcpPi;
        }

        FORCE_INLINE virtual double evaluate_pdf(
            const void*         data,
            const Vector3d&     geometric_normal,
            const Basis3d&      shading_basis,
            const Vector3d&     outgoing,
            const Vector3d&     incoming,
            const int           modes) const
        {
            if (!(modes & Diffuse))
                return 0.0;

            // No reflection below the shading surface.
            const Vector3d& n = shading_basis.get_normal();
            const double cos_in = dot(incoming, n);
            if (cos_in < 0.0)
                return 0.0;

            return cos_in * RcpPi;
        }

      private:
        typedef LommelBRDFInputValues InputValues;
    };

    typedef BSDFWrapper<LommelBRDFImpl> LommelBRDF;
}


//
// LommelBRDFFactory class implementation.
//

const char* LommelBRDFFactory::get_model() const
{
    return Model;
}

const char* LommelBRDFFactory::get_human_readable_model() const
{
    return "Lommel BRDF";
}

DictionaryArray LommelBRDFFactory::get_input_metadata() const
{
    DictionaryArray metadata;

    metadata.push_back(
        Dictionary()
            .insert("name", "reflectance")
            .insert("label", "Reflectance")
            .insert("type", "colormap")
            .insert("entity_types",
                Dictionary()
                    .insert("color", "Colors")
                    .insert("texture_instance", "Textures"))
            .insert("use", "required")
            .insert("default", "0.5"));

    metadata.push_back(
        Dictionary()
            .insert("name", "reflectance_multiplier")
            .insert("label", "Reflectance Multiplier")
            .insert("type", "colormap")
            .insert("entity_types",
                Dictionary().insert("texture_instance", "Textures"))
            .insert("use", "optional")
            .insert("default", "1.0"));

    metadata.push_back(
        Dictionary()
            .insert("name", "thickness")
            .insert("label", "Thickness")
            .insert("type", "colormap")
            .insert("entity_types",
                Dictionary().insert("texture_instance", "Textures"))
            .insert("use", "optional")
            .insert("default", "0.1"));
  

    return metadata;
}

auto_release_ptr<BSDF> LommelBRDFFactory::create(
    const char*         name,
    const ParamArray&   params) const
{
    return auto_release_ptr<BSDF>(new LommelBRDF(name, params));
}

}   // namespace renderer
