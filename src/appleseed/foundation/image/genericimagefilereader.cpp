
//
// This source file is part of appleseed.
// Visit http://appleseedhq.net/ for additional information and resources.
//
// This software is released under the MIT license.
//
// Copyright (c) 2010-2013 Francois Beaune, Jupiter Jazz Limited
// Copyright (c) 2014-2015 Francois Beaune, The appleseedhq Organization
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
#include "genericimagefilereader.h"

// appleseed.foundation headers.
#include "foundation/core/exceptions/exceptionunsupportedfileformat.h"
#include "foundation/image/exrimagefilereader.h"
#include "foundation/image/pngimagefilereader.h"
#include "foundation/utility/string.h"

// Boost headers.
#include "boost/filesystem/path.hpp"

// Standard headers.
#include <string>

using namespace boost;
using namespace std;

namespace foundation
{

//
// GenericImageFileReader class implementation.
//

Image* GenericImageFileReader::read(
    const char*         filename,
    ImageAttributes*    image_attributes)
{
    const filesystem::path filepath(filename);
    const string extension = lower_case(filepath.extension().string());

    if (extension == ".exr")
    {
        EXRImageFileReader reader;
        return reader.read(filename, image_attributes);
    }
    else if (extension == ".png")
    {
        PNGImageFileReader reader;
        return reader.read(filename, image_attributes);
    }
    else
    {
        throw ExceptionUnsupportedFileFormat(filename);
        return 0;
    }
}

}   // namespace foundation
