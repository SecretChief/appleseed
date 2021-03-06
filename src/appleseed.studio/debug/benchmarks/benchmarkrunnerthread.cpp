
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
#include "benchmarkrunnerthread.h"

// appleseed.shared headers.
#include "application/application.h"

// appleseed.foundation headers.
#include "foundation/platform/thread.h"
#include "foundation/utility/autoreleaseptr.h"
#include "foundation/utility/benchmark.h"
#include "foundation/utility/string.h"

// Boost headers.
#include "boost/filesystem/operations.hpp"
#include "boost/filesystem/path.hpp"

// Standard headers.
#include <cassert>
#include <string>

using namespace appleseed::shared;
using namespace boost;
using namespace foundation;
using namespace std;

namespace appleseed {
namespace studio {

//
// BenchmarkRunnerThread class implementation.
//

void BenchmarkRunnerThread::run()
{
    set_current_thread_name("benchmarks");

    auto_release_ptr<XMLFileBenchmarkListener> xmlfile_listener(
        create_xmlfile_benchmark_listener());

    const string xmlfile_name = "benchmark." + get_time_stamp_string() + ".xml";
    const filesystem::path xmlfile_path =
          filesystem::path(Application::get_tests_root_path())
        / "unit benchmarks"
        / "results"
        / xmlfile_name;

    if (!xmlfile_listener->open(xmlfile_path.string().c_str()))
    {
        emit signal_cannot_create_benchmark_file();
        return;
    }

    BenchmarkResult result;
    result.add_listener(xmlfile_listener.get());

    const filesystem::path old_current_path =
        Application::change_current_directory_to_tests_root_path();

    BenchmarkSuiteRepository::instance().run(result);

    filesystem::current_path(old_current_path);

    emit signal_finished();
}

}   // namespace studio
}   // namespace appleseed
