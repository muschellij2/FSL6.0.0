/*  Multimodal Image Segmentation Tool (MIST)  */
/*  Eelke Visser  */

/*  Copyright (c) 2016 University of Oxford  */

/*  Part of FSL - FMRIB's Software Library
    http://www.fmrib.ox.ac.uk/fsl
    fsl@fmrib.ox.ac.uk
    
    Developed at FMRIB (Oxford Centre for Functional Magnetic Resonance
    Imaging of the Brain), Department of Clinical Neurology, Oxford
    University, Oxford, UK
    
    
    LICENCE
    
    FMRIB Software Library, Release 5.0 (c) 2012, The University of
    Oxford (the "Software")
    
    The Software remains the property of the University of Oxford ("the
    University").
    
    The Software is distributed "AS IS" under this Licence solely for
    non-commercial use in the hope that it will be useful, but in order
    that the University as a charitable foundation protects its assets for
    the benefit of its educational and research purposes, the University
    makes clear that no condition is made or to be implied, nor is any
    warranty given or to be implied, as to the accuracy of the Software,
    or that it will be suitable for any particular purpose or for use
    under any specific conditions. Furthermore, the University disclaims
    all responsibility for the use which is made of the Software. It
    further disclaims any liability for the outcomes arising from using
    the Software.
    
    The Licensee agrees to indemnify the University and hold the
    University harmless from and against any and all claims, damages and
    liabilities asserted by third parties (including claims for
    negligence) which arise directly or indirectly from the use of the
    Software or the sale of any products based on the Software.
    
    No part of the Software may be reproduced, modified, transmitted or
    transferred in any form or by any means, electronic or mechanical,
    without the express permission of the University. The permission of
    the University is not required if the said reproduction, modification,
    transmission or transference is done without financial return, the
    conditions of this Licence are imposed upon the receiver of the
    product, and all original and amended source code is included in any
    transmitted product. You may be held legally responsible for any
    copyright infringement that is caused or encouraged by your failure to
    abide by these terms and conditions.
    
    You are not permitted under this Licence to use this Software
    commercially. Use for which any financial return is received shall be
    defined as commercial use, and includes (1) integration of all or part
    of the source code or the Software into a product for sale or license
    by or on behalf of Licensee to third parties or (2) use of the
    Software or any derivative of it for research with the final aim of
    developing software products for sale or license to a third party or
    (3) use of the Software or any derivative of it for research with the
    final aim of developing non-software products for sale or license to a
    third party, or (4) use of the Software to provide any service to an
    external organisation for which payment is received. If you are
    interested in using the Software commercially, please contact Oxford
    University Innovation ("OUI"), the technology transfer company of the
    University, to negotiate a licence. Contact details are:
    Innovation@innovation.ox.ac.uk quoting reference DE/9564. */

#include "plotting.h"
#include <cstdio>
#include <boost/log/trivial.hpp>
#include <iostream>
#include <fstream>
#include <sstream>

extern "C"
{
    #include "gd.h"
    #include "gdc.h"
    #include "gdchart.h"
}

#define WANT_STREAM
#include "newmatio.h"

using namespace NEWMAT;

namespace Plotting
{
    static char ylabelfmt[] = "%f";

    void plotpng(const std::vector<ColumnVector> &data, const::vector<std::size_t> &deltas, int steps,
                 int vertex, const std::string &plotname, const std::string &basename)
    {
        BOOST_LOG_TRIVIAL(debug) << "Creating plot with " << data.size() << " lines and " << data[0].Nrows() << " points ...";

        std::vector<float> floats;
        std::size_t points = data[0].Nrows();

        for (std::size_t j = 0; j < data.size(); j++)
        {
            if (data[j].Nrows() != points)
                throw PlottingException("All lines must have the same length");

            std::size_t delta = 0;
            if (steps)
                delta = deltas[j];

            for (std::size_t i = 0; i < delta; i++)
                floats.push_back(data[j](1));

            for (std::size_t i = 0; i < points; i++)
                floats.push_back(data[j](i + 1));

            for (std::size_t i = 0; i < steps - delta; i++)
                floats.push_back(data[j](points));
        }

        std::string basefilename = basename + "_" + std::to_string(vertex) + "_" + plotname;

        GDC_image_type = GDC_PNG;
        GDC_ylabel_fmt = ylabelfmt;
        FILE *outfile = std::fopen((basefilename + ".png").c_str(), "wb");
        GDC_out_graph(data[0].Nrows() * 15, 400, outfile, GDC_LINE, points + steps, NULL, data.size(), floats.data(), NULL);
        std::fclose(outfile);

        /*
        {
            std::ofstream os(basefilename + ".txt");

            for (std::size_t i = 1; i <= data[0].Nrows(); i++)
                for (std::size_t j = 0; j < data.size(); j++)
                    os << data[j](i) << (j == data.size() - 1 ? "\n" : "\t");
        }
        */

        BOOST_LOG_TRIVIAL(debug) << "... done";
    }

    void plotsqlite(const std::vector<ColumnVector> &data, const::vector<std::size_t> &deltas, int steps,
                    int vertex, const std::string &plotname, const std::string &dbfile, int jobid)
    {
        std::ostringstream yvalstr;

        std::size_t vecs = data.size();
        for (std::size_t j = 0; j < vecs; j++)
        {
            yvalstr << "(";

            std::size_t rows = data[j].Nrows();
            for (std::size_t i = 1; i <= rows; i++)
            {
                yvalstr << data[j](i);

                if (i != rows)
                    yvalstr << ", ";
            }

            yvalstr << ")";

            if (j != vecs - 1)
                yvalstr << ", ";
        }

        std::ostringstream deltastr;
        deltastr << "(";

        for (std::size_t i = 0; i < deltas.size(); i++)
        {
            deltastr << deltas[i];

            if (i != deltas.size() - 1)
                deltastr << ", ";
        }

        deltastr << ")";

        BOOST_LOG_TRIVIAL(debug) << "Opening sqlite database";

        sqlite3 *db;
        int sq3result = sqlite3_open_v2(dbfile.c_str(), &db,
                                                SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE, NULL);

        if (sq3result != SQLITE_OK)
            throw PlottingException("Cannot open sqlite database. Error code: " + std::to_string(sq3result));

        // Wait for up to one hour for the database to be unlocked
        sqlite3_busy_timeout(db, 3600000);

        std::string createtext = std::string("CREATE TABLE IF NOT EXISTS plots "
                               "(jobid INT, vertex INT, name TEXT, yvals TEXT, deltas TEXT, steps INT, "
                               "CONSTRAINT pk PRIMARY KEY (jobid, vertex, name));");

        BOOST_LOG_TRIVIAL(debug) << "Creating table (if non-existent)";

        sqlite3_stmt *createstmt;
        int result = sqlite3_prepare_v2(db, createtext.c_str(), -1, &createstmt, NULL);
        if (result != SQLITE_OK)
            throw PlottingException("sqlite3_prepare_v2() failed for CREATE statement with error code " + std::to_string(result));

        result = sqlite3_step(createstmt);
        if (result != SQLITE_DONE)
            throw PlottingException("sqlite3_step() returned error code " + std::to_string(result)
                                   + " for CREATE statement (SQLITE_DONE was expected)");

        sqlite3_finalize(createstmt);

        BOOST_LOG_TRIVIAL(debug) << "Inserting plot data";

        std::string inserttext = "INSERT INTO plots VALUES (?, ?, ?, ?, ?, ?)";

        sqlite3_stmt *insertstmt;
        result = sqlite3_prepare_v2(db, inserttext.c_str(), -1, &insertstmt, NULL);
        if (result != SQLITE_OK)
            throw PlottingException("sqlite3_prepare_v2() failed for INSERT statement with error code " + std::to_string(result));

        result = sqlite3_bind_int(insertstmt, 1, jobid);
        if (result != SQLITE_OK)
            throw PlottingException("sqlite3_bind_int() failed to bind jobid with error code " + std::to_string(result));

        result = sqlite3_bind_int(insertstmt, 2, vertex);
        if (result != SQLITE_OK)
            throw PlottingException("sqlite3_bind_int() failed to bind vertex with error code " + std::to_string(result));

        result = sqlite3_bind_text(insertstmt, 3, plotname.c_str(), -1, SQLITE_TRANSIENT);
        if (result != SQLITE_OK)
            throw PlottingException("sqlite3_bind_text() failed to bind name with error code " + std::to_string(result));

        result = sqlite3_bind_text(insertstmt, 4, yvalstr.str().c_str(), -1, SQLITE_TRANSIENT);
        if (result != SQLITE_OK)
            throw PlottingException("sqlite3_bind_text() failed to bind yvals with error code " + std::to_string(result));

        result = sqlite3_bind_text(insertstmt, 5, deltastr.str().c_str(), -1, SQLITE_TRANSIENT);
        if (result != SQLITE_OK)
            throw PlottingException("sqlite3_bind_text() failed to bind deltas with error code " + std::to_string(result));

        result = sqlite3_bind_int(insertstmt, 6, steps);
        if (result != SQLITE_OK)
            throw PlottingException("sqlite3_bind_int() failed to bind steps with error code " + std::to_string(result));

        result = sqlite3_step(insertstmt);
        if (result != SQLITE_DONE)
            throw PlottingException("sqlite3_step() returned error code " + std::to_string(result)
                                   + " for INSERT statement (SQLITE_DONE was expected)");

        sqlite3_finalize(insertstmt);

        BOOST_LOG_TRIVIAL(debug) << "Closing database";

        sq3result = sqlite3_close_v2(db);
        if (sq3result != SQLITE_OK)
            throw PlottingException("Cannot close sqlite database. Error code: " + std::to_string(sq3result));
    }
}
