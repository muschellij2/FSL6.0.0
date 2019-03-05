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

#include "profilemixtures.h"
#include "mvnshapemodel.h"
#include "mrfshapemodel.h"
#include "miscmaths/miscmaths.h"
#include "profilepriors.h"
#include "profilefilters.h"
#include "plotting.h"
#include "transformation.h"
#include "newimageall.h"
#include <random>
// #include <boost/date_time/posix_time/posix_time_types.hpp>
#include <boost/log/core.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/utility/setup/console.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>
#include <boost/log/support/date_time.hpp>
#include <boost/serialization/export.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/regex.hpp>
#include <boost/make_shared.hpp>
#include "utils/options.h"
#include <string>
#include <vector>
#include <unordered_map>
#include <cstdlib>
#include "builddate.h"

// TODO: Maybe add normalisation in maskmedians mode (probably doesn't make much difference and adding normalisation
//       modes makes things even more complicated)
// TODO: Look into BOOST_STATIC_ASSERT_UNUSED_ATTRIBUTE warnings (seems like it might be fixed in next boost release)
// TODO FIX: Plot of aligned profiles to png seems to be broken - sqlite is fine
// NB: 'Aligned plot' is shorter than length of profile that is aligned when outputting png!
// NB: There may be issues loading inf/nan/etc from files (but they shouldn't be there in the first place)
// TODO: Show default values in help!

using namespace NEWIMAGE;
using namespace Utilities;

void train(int argc, char *argv[]);
void fit(int argc, char *argv[]);
void dump(int argc, char *argv[]);
void maskmedians(int argc, char *argv[]);

void exiterror(OptionParser &parser, string error)
{
    parser.usage();
    std::cerr << "\nError: " << error << std::endl;

    exit(EXIT_FAILURE);
}

std::vector<std::string> splitarg(std::string arg)
{
    std::vector<std::string> elements;
    boost::algorithm::split(elements, arg, boost::algorithm::is_any_of(","));

    return elements;
}

boost::shared_ptr<ProfilePrior> parsepriorspec(std::string spec)
{
    boost::regex flat("flat\\(([\\d\\.\\-eE]+)\\)");
    boost::regex step("step\\(([\\d\\.\\-eE]+);([\\d\\.\\-eE]+)\\)");
    boost::regex block("block\\(([\\d\\.\\-eE]+);([\\d\\.\\-eE]+);([\\d\\.\\-eE]+)\\)");
    boost::regex ra("ra\\(([\\d\\.\\-eE]+);([\\d\\.\\-eE]+);([\\d\\.\\-eE]+)\\)");
    boost::regex expp("exp\\(([\\d\\.\\-eE]+);([\\d\\.\\-eE]+);([\\d\\.\\-eE]+)\\)");
    boost::regex exp2p("exp2\\(([\\d\\.\\-eE]+);([\\d\\.\\-eE]+);([\\d\\.\\-eE]+);([\\d\\.\\-eE]+)\\)");
    boost::regex doubleexpp("doubleexp\\(([\\d\\.\\-eE]+);([\\d\\.\\-eE]+);([\\d\\.\\-eE]+);([\\d\\.\\-eE]+);([\\d\\.\\-eE]+)\\)");
    boost::regex doubleexp2p("doubleexp2\\(([\\d\\.\\-eE]+);([\\d\\.\\-eE]+);([\\d\\.\\-eE]+);([\\d\\.\\-eE]+);([\\d\\.\\-eE]+);([\\d\\.\\-eE]+)\\)");
    boost::regex sq("sq\\(([\\d\\.\\-eE]+);([\\d\\.\\-eE]+)\\)");

    boost::smatch m;

    if (regex_match(spec, m, flat))
        return boost::make_shared<FlatPrior>(std::stod(m[1]));

    if (regex_match(spec, m, step))
        return boost::make_shared<SimpleEdgePrior>(std::stod(m[1]), std::stod(m[2]));

    if (regex_match(spec, m, block))
        return boost::make_shared<BlockPrior>(std::stod(m[1]), std::stod(m[2]), std::stod(m[3]));

    if (regex_match(spec, m, ra))
        return boost::make_shared<RightAngledPrior>(std::stod(m[1]), std::stod(m[2]), std::stod(m[3]));

    if (regex_match(spec, m, expp))
        return boost::make_shared<ExponentialPrior>(std::stod(m[1]), std::stod(m[2]), std::stod(m[3]));

    if (regex_match(spec, m, exp2p))
        return boost::make_shared<Exponential2Prior>(std::stod(m[1]), std::stod(m[2]), std::stod(m[3]),
                std::stod(m[4]));
    
    if (regex_match(spec, m, doubleexpp))
        return boost::make_shared<DoubleExponentialPrior>(std::stod(m[1]), std::stod(m[2]), std::stod(m[3]),
                std::stod(m[4]), std::stod(m[5]));
    
    if (regex_match(spec, m, doubleexp2p))
        return boost::make_shared<DoubleExponential2Prior>(std::stod(m[1]), std::stod(m[2]), std::stod(m[3]),
                std::stod(m[4]), std::stod(m[5]), std::stod(m[6]));

    if (regex_match(spec, m, sq))
        return boost::make_shared<ParabolicPrior>(std::stod(m[1]), std::stod(m[2]));

    return nullptr;
}

int main(int argc, char *argv[])
{
    boost::log::add_common_attributes();
    boost::log::add_console_log(std::cerr, boost::log::keywords::format = boost::log::expressions::stream
            << "[" << boost::log::expressions::format_date_time<boost::posix_time::ptime>("TimeStamp", "%Y-%m-%d %H:%M:%S.%f") << "]"
            << " [" << boost::log::trivial::severity << "]"
            << "\t" << boost::log::expressions::smessage);

    BOOST_LOG_TRIVIAL(info) << "Multimodal image segmentation tool (MIST), built " << BUILDDATE ".";

    try
    {
        OptionParser parser("Multimodal image segmentation tool (MIST)\nCopyright(c) 2015, University of Oxford (Eelke Visser)", "mist [options] <train|fit|maskmedians> <options>");

        Option<bool> globalhelp("-h,--help", false, "Show help", false, no_argument);
        parser.add(globalhelp);

        Option<string> loglevel("--loglevel", "info", "Set logging level", false, requires_argument);
        parser.add(loglevel);

        int nonopt;
        try
        {
            nonopt = parser.parse_command_line(argc, argv, 0, true);
        }
        catch (X_OptionError &e)
        {
            exiterror(parser, e.what());
        }

        if (loglevel.set())
        {
            int level;

            if (loglevel.value() == "trace")
                level = boost::log::trivial::trace;
            else if (loglevel.value() == "debug")
                level = boost::log::trivial::debug;
            else if (loglevel.value() == "info")
                level = boost::log::trivial::info;
            else if (loglevel.value() == "warning")
                level = boost::log::trivial::warning;
            else if (loglevel.value() == "error")
                level = boost::log::trivial::error;
            else if (loglevel.value() == "fatal")
                level = boost::log::trivial::fatal;
            else
                exiterror(parser, "Invalid log level");

            boost::log::core::get()->set_filter(boost::log::trivial::severity >= level);
        }
        else
            boost::log::core::get()->set_filter(boost::log::trivial::severity >= boost::log::trivial::info);

        if (globalhelp.value() == true)
            parser.usage();
        else if (nonopt == argc)
            exiterror(parser, "No mode specified (type 'mist <mode>' to get help on using that mode).");
        else
        {
            // Include mode verb -> same counting convention as global args
            string verb(argv[nonopt]);
            argc -= nonopt;
            argv += nonopt;

            if (verb == "train")
                train(argc, argv);
            else if (verb == "fit")
                fit(argc, argv);
            else if (verb == "dump")
                dump(argc, argv);
            else if (verb == "maskmedians")
                maskmedians(argc, argv);
            else
                exiterror(parser, std::string("No such mode: ") + verb);
        }
    }
    catch (const Exception &e)
    {
        BOOST_LOG_TRIVIAL(fatal) << "Caught exception: " << e.what();

        return EXIT_FAILURE;
    }
    catch (const std::exception &e)
    {
        BOOST_LOG_TRIVIAL(fatal) << "Caught exception: " << e.what();

        return EXIT_FAILURE;
    }

    BOOST_LOG_TRIVIAL(info) << "Done";

    return EXIT_SUCCESS;
}

void train(int argc, char *argv[])
{
    OptionParser parser("Multimodal subcortical segmentation: Training mode", "mist [options] train <options> <subjects ...>");

    Option<string> shapefile("--shape", "", "Filename of reference shape", true, requires_argument);
    parser.add(shapefile);

    Option<string> modalitynames("--modalitynames", "", "Comma-separated list of modality names", true, requires_argument);
    parser.add(modalitynames);

    Option<string> modalityfiles("--modalityimages", "", "Comma-separated list of images relative to subject path", true, requires_argument);
    parser.add(modalityfiles);

    Option<string> affinefile("--affine", "", "Filename of affine transform relative to subject path", false, requires_argument);
    parser.add(affinefile);

    Option<string> warpfile("--warp", "", "Filename of warp field coefficients relative to subject path", false, requires_argument);
    parser.add(warpfile);

    Option<string> vertices("--vertices", "", "Fit specified vertices (comma-separated list; use --out to specify base filename)", false, requires_argument);
    parser.add(vertices);

    Option<string> loadvertices("--loadvertices", "", "Do not fit vertex models but load them from files; specify basename", false, requires_argument);
    parser.add(loadvertices);

    Option<string> outfile("--out", "", "Output filename (serialization archive - either full model or base filename for single vertex models)", true, requires_argument);
    parser.add(outfile);

    Option<int> profilepoints("--profilepoints", 10, "Number of points to sample along normal at each vertex", false, requires_argument);
    parser.add(profilepoints);

    Option<float> profilespacing("--profilespacing", 0.5, "Spacing between these points in mm", false, requires_argument);
    parser.add(profilespacing);

    Option<string> smoothness("--smoothness", "", "Sigma of smoothing kernel for covariance matrix (one value per modality, comma-separated)",
                              true, requires_argument);
    parser.add(smoothness);

    Option<int> components("--components", 1 , "Number of componenents in mixture models", false, requires_argument);
    parser.add(components);

    Option<int> steps("--steps", 20, "Number of steps (+1) from minimum to maximum displacement", false, requires_argument);
    parser.add(steps);

    Option<string> priormeans("--priormeans", "", "Comma-separated list of prior mean specifications. Valid specifications are: flat(x), step(x;y), block(w;x;y), ra(w;x;y), exp(t;x;y), exp2(t;x;y;z), doubleexp(t1;t2;x;y;z), doubleexp2(t1;t2;x;y;z;w), sq(x;y)",
                              true, requires_argument);
    parser.add(priormeans);

    Option<string> priorcovcoefs("--priorcovcoefs", "", "Comma-separated list of prior covariance coefficient specifications (beta is resulting covariance matrix multiplied by alpha). Valid specifications are: flat(x), step(x;y), block(w;x;y), ra(w;x;y), exp(t;x;y), exp2(t;x;y;z), doubleexp(t1;t2;x;y;z), doubleexp2(t1;t2;x;y;z;w), sq(x;y)",
                              true, requires_argument);
    parser.add(priorcovcoefs);

    Option<int> profilen0("--profilen0", 10, "N0 for profile model", false, requires_argument);
    parser.add(profilen0);

    Option<float> profilealpha("--profilealpha", -3.0, "Alpha for profile model (< 0 to specify wrt (dims-1)/2)", false, requires_argument);
    parser.add(profilealpha);

    Option<float> mixingalpha("--mixingalpha", 2.0, "Alpha for dirichlet prior (all components)", false, requires_argument);
    parser.add(mixingalpha);

    Option<float> deltastdev("--deltastdev", 5.0, "Standard deviation of delta in steps (float)", false, requires_argument);
    parser.add(deltastdev);

    Option<float> ftol("--ftol", 1e-10, "Fractional tolerance on profile model log probability", false, requires_argument);
    parser.add(ftol);

    Option<int> ftolignore("--ftolignore", 4, "Number of times to ignore reaching fractional tolerance", false, requires_argument);
    parser.add(ftolignore);

    Option<int> maxeval("--maxeval", 1000, "Maximum number of evaluations of the above", false, requires_argument);
    parser.add(maxeval);

    Option<int> maxretry("--maxretry", 0, "Maximum number of retries for optimisation", false, requires_argument);
    parser.add(maxretry);

    HiddenOption<bool> usemcmc("--usemcmc", false, "Generate MCMC samples for profile models", false, no_argument);
    parser.add(usemcmc);

    Option<string> normalisation("--normalisation", "", "Set normalisation modes for all modalities (modes are none, add, mul; default is mul)", false, requires_argument);
    parser.add(normalisation);

    Option<string> normexclusion("--normexclusion", "", "Exclusion mask for intensity normalisation", false, requires_argument);
    parser.add(normexclusion);

    HiddenOption<string> filter("--filter", "", "Filter for each modality (none, deriv, deriv+, or deriv-)", false, requires_argument);
    parser.add(filter);

    Option<bool> usemrf("--usemrf", false, "Use MRF shape model instead of the MVN one", false, no_argument);
    parser.add(usemrf);

    Option<float> mrfweight("--mrfweight", 1.0, "Weight parameter for MRF", false, requires_argument);
    parser.add(mrfweight);

    Option<float> mrfmeanfrac("--mrfmeanfrac", 0.0, "Relative (fractional) weight of mean term in MRF", false, requires_argument);
    parser.add(mrfmeanfrac);

    Option<int> shapen0("--shapen0", 10, "Prior N0 for shape model", false, requires_argument);
    parser.add(shapen0);

    Option<float> shapealpha("--shapealpha", -3.0, "Alpha for shape model (< 0 to specify wrt (dims-1)/2)", false, requires_argument);
    parser.add(shapealpha);

    Option<float> shapecovstdev("--shapecovstdev", 0.05, "Amplitude parameter for displacement covariance matrix", false, requires_argument);
    parser.add(shapecovstdev);

    Option<float> shapecovwidth("--shapecovwidth", 0.6, "Width parameter for displacement covariance matrix", false, requires_argument);
    parser.add(shapecovwidth);

    Option<string> plotprefix("--plotprefix", "", "Write plots; specify prefix", false, requires_argument);
    parser.add(plotprefix);

    Option<string> plotsqlite("--plotsqlite", "", "Output plot data to table plots (id, vert, data). Specify filename and unique ID. NB: Database location should support locking!", false, requires_argument);
    parser.add(plotsqlite);

    Option<bool> verbhelp("-h,--help", false, "Show help for this mode", false, no_argument);
    parser.add(verbhelp);

    int nextarg;
    try
    {
        nextarg = parser.parse_command_line(argc, argv, 0, true);
    }
    catch (X_OptionError &e)
    {
        exiterror(parser, e.what());
    }

    if (verbhelp.set())
    {
        parser.usage();
        return;
    }

    if (!parser.check_compulsory_arguments(true))
        exiterror(parser, "Not all compulsory arguments were set");

    if (nextarg == argc)
        exiterror(parser, "No subjects specified");

    if (affinefile.set() == warpfile.set())
        exiterror(parser, "Either --affine or --warp should be set");

    if (profilepoints.value() % 2 || steps.value() % 2)
        exiterror(parser, "Both --profilepoints and --steps should be even");

    std::vector<std::string> names = splitarg(modalitynames.value());
    std::vector<std::string> files = splitarg(modalityfiles.value());
    std::vector<std::string> sigmas = splitarg(smoothness.value());
    std::vector<std::string> priormeanspecs = splitarg(priormeans.value());
    std::vector<std::string> priorcovcoefspecs = splitarg(priorcovcoefs.value());

    if (files.size() != names.size()
            || sigmas.size() != names.size())
        exiterror(parser, "The number of modalities should be the same for all argument that specify a list");

    if (priormeanspecs.size() != names.size() * components.value()
            || priorcovcoefspecs.size() != priormeanspecs.size())
        exiterror(parser, "Mean and covariance coefficient priors should be specified for each component in each modality (grouped by modality)");

    if (!usemrf.set() && (mrfweight.set() || mrfmeanfrac.set()))
        exiterror(parser, "MRF parameters are meaningless without --usemrf");

    if (usemrf.set() && (shapecovwidth.set() || shapecovstdev.set() || shapealpha.set() || shapen0.set()))
        exiterror(parser, "MVN parameters are meaningless when using --usemrf");

    std::unordered_map<std::string, std::vector<std::string> > trainingdata;
    std::vector<boost::shared_ptr<const Transformation> > transformations;
    std::vector<std::string> normalisationmasks;

    BOOST_LOG_TRIVIAL(info) << "Loading transformations";
    for (; nextarg != argc; nextarg++)
    {
        std::string subj(argv[nextarg]);

        for (auto n = names.cbegin(), f = files.cbegin(); n != names.cend(); n++, f++)
            trainingdata[*n].push_back(subj + "/" + *f);

        if (affinefile.set())
            transformations.push_back(boost::make_shared<AffineTransformation>(subj + "/" + affinefile.value()));
        else
            transformations.push_back(boost::make_shared<NonlinearTransformation>(subj + "/" + warpfile.value(), true));

        if (normexclusion.set())
            normalisationmasks.push_back(subj + "/" + normexclusion.value());
    }

    // NB: This can fail with 'No file specified' when wrong versions of the VTK libs are loaded
    Shape shape(shapefile.value(), "reference");

    boost::shared_ptr<ProfileMixtures> pm = boost::make_shared<ProfileMixtures>(
                        names, profilepoints.value() + steps.value(), profilepoints.value(), components.value());

    pm->SetN0(profilen0.value());

    if (profilealpha.value() < 0.0)
        pm->SetAlpha((profilepoints.value() + steps.value() - 1.0) / 2.0 - profilealpha.value());
    else
        pm->SetAlpha(profilealpha.value());

    pm->SetDeltaStdev(deltastdev.value());
    pm->SetMixingAlpha(std::vector<double>(components.value(), mixingalpha.value()));

    for (std::size_t m = 0; m < names.size(); m++)
    {
        for (int c = 0; c < components.value(); c++)
        {
            int ind = m * components.value() + c;

            boost::shared_ptr<ProfilePrior> mp = parsepriorspec(priormeanspecs[ind]);
            if (!mp)
                exiterror(parser, "Invalid prior specified for mean");

            pm->SetComponentPriorMean(names[m], c, mp, profilespacing.value());

            boost::shared_ptr<ProfilePrior> ccp = parsepriorspec(priorcovcoefspecs[ind]);
            if (!ccp)
                exiterror(parser, "Invalid prior specified for covariance coefficients");

            pm->SetComponentPriorCovarianceCoefs(names[m], c, ccp, profilespacing.value());
        }

        pm->SetSmoothness(names[m], std::stod(sigmas[m]));
    }

    pm->SetFTolerance(ftol.value());
    pm->SetFToleranceIgnoreCount(ftolignore.value());
    pm->SetMaxEvaluations(maxeval.value());
    pm->SetMaxRetries(maxretry.value());
    pm->SetUseMCMC(usemcmc.value());

    boost::shared_ptr<ShapeModel> model = nullptr;

    if (usemrf.set())
    {
        auto sm = boost::make_shared<MRFShapeModel>(
                        shape, names, profilepoints.value(), profilespacing.value(), normexclusion.set());

        // Don't think this is used for training; number of training iterations is hardcoded for now
        // (and it is reset below for the final segmentation)
        sm->SetMaxIterations(100);

        sm->SetWeight(mrfweight.value());
        sm->SetMeanFraction(mrfmeanfrac.value());

        model = sm;
    }
    else
    {
        auto sm = boost::make_shared<MVNShapeModel>(
                        shape, names, profilepoints.value(), profilespacing.value(), normexclusion.set());

        sm->UseShapeVariability(true);
        sm->SetShapeCovariancePrior(shapecovstdev.value(), shapecovwidth.value());
        sm->SetShapeN0(shapen0.value());
        if (shapealpha.value() < 0.0)
            sm->SetShapeAlpha((shape.GetNumberOfVertices() - 1.0) / 2.0 - shapealpha.value());
        else
            sm->SetShapeAlpha(shapealpha.value());

        model = sm;
    }

    if (normalisation.set())
    {
        std::vector<std::string> normmodes = splitarg(normalisation.value());

        if (normmodes.size() != names.size())
            exiterror(parser, "Normalisation mode must be specified for all modalities if using --normalisation");

        for (int m = 0; m < names.size(); m++)
        {
            if (normmodes[m] == "none")
                model->SetNormalise(names[m], ShapeModel::NormalisationMode::None);
            else if (normmodes[m] == "add")
                model->SetNormalise(names[m], ShapeModel::NormalisationMode::Additive);
            else if (normmodes[m] == "mul")
                model->SetNormalise(names[m], ShapeModel::NormalisationMode::Multiplicative);
            else
                exiterror(parser, std::string("Invalid normalisation mode: ") + normmodes[m]);
        }
    }
    else
        for (int m = 0; m < names.size(); m++)
            model->SetNormalise(names[m], ShapeModel::NormalisationMode::Multiplicative);

    if (filter.set())
    {
        std::vector<std::string> filts = splitarg(filter.value());

        if (filts.size() != names.size())
            exiterror(parser, "A filter (or 'none') should be specified for each modality when using the --filter option");

        for (auto nm = names.cbegin(), filt = filts.cbegin(); nm != names.cend(); nm++, filt++)
        {
            if (*filt == "none")
                model->SetFilter(*nm, boost::make_shared<IdentityFilter>());
            else if (*filt == "deriv")
                model->SetFilter(*nm, boost::make_shared<DerivativeFilter>(DerivativeFilter::Full));
            else if (*filt == "deriv+")
                model->SetFilter(*nm, boost::make_shared<DerivativeFilter>(DerivativeFilter::Positive));
            else if (*filt == "deriv-")
                model->SetFilter(*nm, boost::make_shared<DerivativeFilter>(DerivativeFilter::Negative));
            else
                exiterror(parser, "Unknown filter specified");
        }
    }

    model->CreateVertexModels(pm);

    Plotting::plotfunc pfunc = nullptr;

    if (plotprefix.set() && plotsqlite.set())
        exiterror(parser, "--plotprefix and --plotsqlite are mutually exclusive");

    if (plotprefix.set())
        pfunc = std::bind(&Plotting::plotpng,
                          std::placeholders::_1,
                          std::placeholders::_2,
                          std::placeholders::_3,
                          std::placeholders::_4,
                          std::placeholders::_5,
                          plotprefix.value());

    if (plotsqlite.set())
    {
        char *jobidstr = std::getenv("JOB_ID");
        int jobid = jobidstr == NULL ? 0 : std::atoi(jobidstr);

        pfunc = std::bind(&Plotting::plotsqlite,
                          std::placeholders::_1,
                          std::placeholders::_2,
                          std::placeholders::_3,
                          std::placeholders::_4,
                          std::placeholders::_5,
                          plotsqlite.value(),
                          jobid);
    }

    if (vertices.set())
    {
        std::vector<int> verts;
        for (string &s : splitarg(vertices.value()))
            verts.push_back(std::atoi(s.c_str()));

        std::vector<boost::shared_ptr<ProfileModel> > pms = model->TrainVertices(trainingdata, transformations,
                                                                                 normalisationmasks, verts, pfunc);

        for (int i = 0; i < verts.size(); i++)
        {
            std::ofstream os(outfile.value() + "_vertex" + std::to_string(verts[i]) + ".txt");
            boost::archive::text_oarchive archive(os);
            archive << pms[i];
        }
    }
    else
    {
        if (loadvertices.set())
        {
            std::vector<std::string> vertexmodels;
            for (int i = 0; i < model->GetShape()->GetNumberOfVertices(); i++)
                vertexmodels.push_back(loadvertices.value() + "_vertex" + std::to_string(i) + ".txt");

            model->LoadVertexModelsAndTrain(vertexmodels, trainingdata, transformations, normalisationmasks);
        }
        else
            model->Train(trainingdata, transformations, normalisationmasks);

        std::ofstream os(outfile.value());
        boost::archive::text_oarchive archive(os);

        // Serialize through base class pointer to get type info into archive
        boost::shared_ptr<ShapeModel> baseptr = model;
        archive << baseptr;

        if (plotprefix.set())
            model->WritePlots(pfunc);
    }
}

void fit(int argc, char *argv[])
{
    OptionParser parser("Multimodal subcortical segmentation: Fitting mode", "mist [options] fit <options>");

    Option<string> modelfile("--model", "", "Model filename", true, requires_argument);
    parser.add(modelfile);

    Option<string> modalitynames("--modalitynames", "", "List of modality names", true, requires_argument);
    parser.add(modalitynames);

    Option<string> modalityfiles("--modalityimages", "", "List of images", true, requires_argument);
    parser.add(modalityfiles);

    Option<string> affinefile("--affine", "", "Filename of affine transform", false, requires_argument);
    parser.add(affinefile);

    Option<string> warpfile("--warp", "", "Filename of warp field coefficients relative to subject path", false, requires_argument);
    parser.add(warpfile);

    Option<string> normexclusion("--normexclusion", "", "Exclusion mask for intensity normalisation", false, requires_argument);
    parser.add(normexclusion);

    Option<float> cilevel("--cilevel", 0.9, "Credible interval level", false, requires_argument);
    parser.add(cilevel);

    Option<int> iterations("--iterations", 750, "Number of Gibbs iterations", false, requires_argument);
    parser.add(iterations);

    Option<bool> usegibbs("--usegibbs", false, "Use Gibbs sampling for MVN model", false, no_argument);
    parser.add(usegibbs);

    Option<int> burnin("--burnin", 250, "Number of burn-in iterations for Gibbs", false, requires_argument);
    parser.add(burnin);

    HiddenOption<float> resetmrfweight("--resetmrfweight", 0.0, "Reset the MRF weight parameter (i.e. ignore the value in the model file)", false, requires_argument);
    parser.add(resetmrfweight);

    HiddenOption<float> resetmrfmeanfrac("--resetmrfmeanfrac", 0.0, "Reset relative weight of mean for MRF (i.e. ignore the value in the model file)", false, requires_argument);
    parser.add(resetmrfmeanfrac);

    HiddenOption<bool> disablemcmc("--disablemcmc", false, "Disable use of MCMC samples in profile model even in the model contains them", false, no_argument);
    parser.add(disablemcmc);

    Option<bool> nointersectremove("--nointersectremove", false, "Do not remove self-intersections", false, no_argument);
    parser.add(nointersectremove);

    HiddenOption<bool> nofit("--nofit", false, "Do not fit model; just return the undeformed mesh", false, no_argument);
    parser.add(nofit);

    Option<string> outbase("--outbase", "", "Output base filename for displacements", true, requires_argument);
    parser.add(outbase);

    Option<string> outregmat("--outregmat", "", "Affine matrix for registered output (this is the matrix that registers the desired output modality to the working coordinate grid)", false, requires_argument);
    parser.add(outregmat);

    Option<string> outregref("--outregref", "", "Reference image for registered output", false, requires_argument);
    parser.add(outregref);

    Option<bool> voxelvertices("--voxelvertices", false, "Use voxel vertices instead of centers to test whether inside mesh (produces larger mesh)", false, no_argument);
    parser.add(voxelvertices);

    Option<int> setseed("--setseed", 0, "Set random seed. The default is to use a fixed seed of 0", false, requires_argument);
    parser.add(setseed);

    Option<bool> seedrandomly("--seedrandomly", false, "Use a random seed. The default is to use a fixed seed of 0", false, no_argument);
    parser.add(seedrandomly);

    Option<bool> verbhelp("-h,--help", false, "Show help for this mode", false, no_argument);
    parser.add(verbhelp);

    int nextarg;
    try
    {
        nextarg = parser.parse_command_line(argc, argv, 0, true);
    }
    catch (X_OptionError &e)
    {
        exiterror(parser, e.what());
    }

    if (verbhelp.set())
    {
        parser.usage();
        return;
    }

    if (!parser.check_compulsory_arguments(true))
        exiterror(parser, "Not all compulsory arguments were set");

    if (nextarg != argc)
        exiterror(parser, "No non-option arguments expected");

    if (affinefile.set() == warpfile.set())
        exiterror(parser, "Either --affine or --warp should be set");

    if (cilevel.value() < 0.0 || cilevel.value() > 1.0)
        exiterror(parser, "Credible interval level should be between 0.0 and 1.0");

    if (outregmat.set() != outregref.set())
        exiterror(parser, "--outregmat and --outregref need to be specified together.");

    if (setseed.set() && seedrandomly.set())
        exiterror(parser, "--setseed and --seedrandomly are mutually exclusive.");

    boost::shared_ptr<const Transformation> transformation;
    if (affinefile.set())
        transformation = boost::make_shared<const AffineTransformation>(affinefile.value());
    else
        transformation = boost::make_shared<const NonlinearTransformation>(warpfile.value(), true);

    std::vector<std::string> names = splitarg(modalitynames.value());
    std::vector<std::string> files = splitarg(modalityfiles.value());

    std::unordered_map<std::string, std::string> data;
    for (auto n = names.cbegin(), f = files.cbegin(); n != names.cend(); n++, f++)
        data[*n] = *f;

    BOOST_LOG_TRIVIAL(info) << "Loading model";
    boost::shared_ptr<ShapeModel> model = nullptr;
    {
        std::ifstream is(modelfile.value());
        boost::archive::text_iarchive archive(is);
        archive >> model;
    }
    BOOST_LOG_TRIVIAL(info) << "Model load complete";

    int verts = model->GetShape()->GetNumberOfVertices();

    model->SetMaxIterations(100);

    if (auto gibbsmodel = boost::dynamic_pointer_cast<GibbsShapeModel>(model))
    {
        BOOST_LOG_TRIVIAL(info) << "Model type is Gibbs";

        if (resetmrfweight.set() || resetmrfmeanfrac.set())
            exiterror(parser, "Cannot reset MRF parameters as this is not an MRF model");

        gibbsmodel->SetGibbsIterations(iterations.value());
        gibbsmodel->SetGibbsBurnIn(burnin.value());

        gibbsmodel->UseGibbs(usegibbs.value());

        if (setseed.set())
            gibbsmodel->Reseed(setseed.value());
        else if (!seedrandomly.set())
            gibbsmodel->Reseed(0);
    }
    else if (auto mrfmodel = boost::dynamic_pointer_cast<MRFShapeModel>(model))
    {
        BOOST_LOG_TRIVIAL(info) << "Model type is MRF";

        if (usegibbs.set() || setseed.set() || seedrandomly.set() || iterations.set() || burnin.set())
            exiterror(parser, "Cannot set Gibbs parameters as this is not a Gibbs model");

        if (resetmrfweight.set())
            mrfmodel->SetWeight(resetmrfweight.value());

        if (resetmrfmeanfrac.set())
            mrfmodel->SetMeanFraction(resetmrfmeanfrac.value());
    }
    else
        BOOST_LOG_TRIVIAL(info) << "Model type is unknown";

    if (disablemcmc.set())
    {
        for (int i = 0; i < verts; i++)
        {
            if (auto pm = boost::dynamic_pointer_cast<ProfileMixtures>(model->GetVertexModel(i)))
                pm->SetUseMCMC(false);
            else
                exiterror(parser, "Could not downcast profile model to ProfileMixtures");
        }
    }

    ColumnVector displacements(verts), cilower(verts), ciupper(verts);
    if (nofit.set())
    {
        BOOST_LOG_TRIVIAL(warning) << "Not fitting to data (disabled on command line)";

        displacements = 0.0;
        cilower = 0.0;
        ciupper = 0.0;
    }
    else
        displacements = model->Fit(data, transformation, normexclusion.value(), cilower, ciupper, cilevel.value());

    Shape finalshape(*model->GetShape(), transformation);
    for (int i = 0; i < finalshape.GetNumberOfVertices(); i++)
        finalshape.SetDisplacement(i, displacements(i + 1));

    if (!nointersectremove.set())
    {
        Shape newshape = finalshape.RemoveSelfIntersections(0.1);

        int changed = 0;
        for (int vert = 0; vert < newshape.GetNumberOfVertices(); vert++)
        {
            if (displacements(vert + 1) != newshape.GetDisplacement(vert))
                changed++;

            displacements(vert + 1) = newshape.GetDisplacement(vert);
        }

        BOOST_LOG_TRIVIAL(info) << "Self-intersection removal adjusted " << changed << " vertices";

        finalshape = newshape;
    }

    BOOST_LOG_TRIVIAL(info) << "Writing primary output files";

    finalshape.WritePolyData(outbase.value() + "_shape.mim");

    {
        ofstream dispout(outbase.value() + "_displacements.txt");
        for (std::size_t i = 1; i <= displacements.Nrows(); i++)
            dispout << displacements(i) << std::endl;
    }

    if (usegibbs.set())
    {
        ofstream cilout(outbase.value() + "_cilower.txt");
        for (std::size_t i = 1; i <= cilower.Nrows(); i++)
            cilout << cilower(i) << std::endl;

        ofstream ciuout(outbase.value() + "_ciupper.txt");
        for (std::size_t i = 1; i <= ciupper.Nrows(); i++)
            ciuout << ciupper(i) << std::endl;
    }

    BOOST_LOG_TRIVIAL(info) << "Writing native mask";

    volume<double> examplevol;
    read_volume(examplevol, data.cbegin()->second);

    volume<char> outvol(examplevol.xsize(), examplevol.ysize(), examplevol.zsize());
    copybasicproperties(examplevol, outvol);

    finalshape.ToVolume(outvol, voxelvertices.set());

    write_volume(outvol, outbase.value() + "_mask");

    if (outregmat.set())
    {
        BOOST_LOG_TRIVIAL(info) << "Writing non-native masks and shape";

        volume<double> refvol;
        read_volume(refvol, outregref.value());
        
        volume<char> outvol(refvol.xsize(), refvol.ysize(), refvol.zsize());
        copybasicproperties(refvol, outvol);

        // TODO: Check for nullptr
        boost::shared_ptr<const Transformation> regtransformation;

        if (affinefile.set())
            regtransformation = boost::make_shared<const AffineTransformation>(
                        AffineTransformation(outregmat.value()).Concatenate(
                        *boost::dynamic_pointer_cast<const AffineTransformation>(transformation)));
        else
            regtransformation = boost::make_shared<const NonlinearTransformation>(
                        AffineTransformation(outregmat.value()).Concatenate(
                        *boost::dynamic_pointer_cast<const NonlinearTransformation>(transformation)));

        Shape regshape(*model->GetShape(), regtransformation);
        for (int i = 0; i < regshape.GetNumberOfVertices(); i++)
            regshape.SetDisplacement(i, displacements(i + 1));

        regshape.WritePolyData(outbase.value() + "_shape_reg.mim");

        regshape.ToVolume(outvol, voxelvertices.set());
        write_volume(outvol, outbase.value() + "_mask_reg");
    }
}

void writecolumnvector(const ColumnVector &vector, const std::string &filename)
{
    std::ofstream os(filename);

    for (std::size_t i = 1; i <= vector.Nrows(); i++)
        os << vector(i) << "\n";
}

void dump(int argc, char *argv[])
{
    OptionParser parser("Multimodal subcortical segmentation: Model dump", "mist [options] dump <options>");

    Option<string> modelfile("--model", "", "Model filename", true, requires_argument);
    parser.add(modelfile);

    Option<string> outroot("--outroot", "", "Output basename", true, requires_argument);
    parser.add(outroot);

    Option<bool> verbhelp("-h,--help", false, "Show help for this mode", false, no_argument);
    parser.add(verbhelp);

    int nextarg;
    try
    {
        nextarg = parser.parse_command_line(argc, argv, 0, true);
    }
    catch (X_OptionError &e)
    {
        exiterror(parser, e.what());
    }

    if (verbhelp.set())
    {
        parser.usage();
        return;
    }

    if (!parser.check_compulsory_arguments(true))
        exiterror(parser, "Not all compulsory arguments were set");

    if (nextarg != argc)
        exiterror(parser, "No non-option arguments expected");

    BOOST_LOG_TRIVIAL(info) << "Loading model";
    boost::shared_ptr<ShapeModel> model;
    {
        std::ifstream is(modelfile.value());
        boost::archive::text_iarchive archive(is);
        archive >> model;
    }
    BOOST_LOG_TRIVIAL(info) << "Model load complete";

    {
        std::ofstream os(outroot.value() + "_profilepoints.txt");
        os << model->GetProfilePoints() << "\n";
    }

    {
        std::ofstream os(outroot.value() + "_profilespacing.txt");
        os << model->GetProfileSpacing() << "\n";
    }

    for (int v = 0; v < model->GetShape()->GetNumberOfVertices(); v++)
    {
        // TODO: Handle bad cast
        auto vertexmodel = boost::dynamic_pointer_cast<const ProfileMixtures>(model->GetVertexModel(v));

        auto priormeans = vertexmodel->GetComponentPriorMeans();
        auto priorcovcoefs = vertexmodel->GetComponentPriorCovarianceCoefs();
        auto means = vertexmodel->GetComponentMeans();
        auto covcoefs = vertexmodel->GetComponentCovarianceCoefs();
        auto mixingcoefs = vertexmodel->GetFullMixingCoefs();
        auto smoothness = vertexmodel->GetSmoothness();

        for (const auto &mn : model->GetModalityNames())
        {
            {
                std::ofstream os(outroot.value() + "_" + std::to_string(v) + "_" + mn + "_mixingcoefs.txt");
                for (double val : mixingcoefs[mn])
                    os << val << "\n";
            }

            {
                std::ofstream os(outroot.value() + "_" + std::to_string(v) + "_" + mn + "_smoothness.txt");
                os << smoothness[mn] << "\n";
            }

            for (int c = 0; c < vertexmodel->GetNumberOfComponents(); c++)
            {
                writecolumnvector(priormeans[mn][c],
                                  outroot.value() + "_" + std::to_string(v) + "_" + mn
                                 + "_" + std::to_string(c) + "_priormean.txt");
                writecolumnvector(priorcovcoefs[mn][c],
                                  outroot.value() + "_" + std::to_string(v) + "_" + mn
                                  + "_" + std::to_string(c) + "_priorcovcoefs.txt");
                writecolumnvector(means[mn][c],
                                  outroot.value() + "_" + std::to_string(v) + "_" + mn
                                  + "_" + std::to_string(c) + "_mean.txt");
                writecolumnvector(covcoefs[mn][c],
                                  outroot.value() + "_" + std::to_string(v) + "_" + mn
                                  + "_" + std::to_string(c) + "_covcoefs.txt");

                Matrix G = ProfileMixtureProbability::GetSmoothingMatrix(smoothness[mn] * smoothness[mn], priorcovcoefs[mn][c].Nrows());

                Matrix Cprior = G * priorcovcoefs[mn][c].AsDiagonal() * G;
                ColumnVector priorstdev(priorcovcoefs[mn][c].Nrows());

                for (std::size_t i = 1; i <= priorcovcoefs[mn][c].Nrows(); i++)
                    priorstdev(i) = std::sqrt(Cprior(i, i));

                Matrix C = G * covcoefs[mn][c].AsDiagonal() * G;
                ColumnVector stdev(covcoefs[mn][c].Nrows());

                for (std::size_t i = 1; i <= covcoefs[mn][c].Nrows(); i++)
                    stdev(i) = std::sqrt(C(i, i));

                writecolumnvector(priorstdev,
                                  outroot.value() + "_" + std::to_string(v) + "_" + mn
                                  + "_" + std::to_string(c) + "_priorstdev.txt");
                writecolumnvector(stdev,
                                  outroot.value() + "_" + std::to_string(v) + "_" + mn
                                  + "_" + std::to_string(c) + "_stdev.txt");
            }
        }
    }
}

void maskmedians(int argc, char *argv[])
{
    OptionParser parser("Multimodal subcortical segmentation: Mask medians", "mist [options] maskmedians <options> <subjects ...>");

    Option<string> imagefiles("--images", "", "Comma-separated list of relative image filenames", true, requires_argument);
    parser.add(imagefiles);

    Option<string> maskfiles("--masks", "", "Comma-separated list of masks", false, requires_argument);
    parser.add(maskfiles);

    Option<string> nativemaskfiles("--nativemasks", "", "Comma-separated list of native masks (relative filenames, medians are output after standard space masks)", false, requires_argument);
    parser.add(nativemaskfiles);

    Option<string> affinefile("--affine", "", "Relative filename of affine matrix", false, requires_argument);
    parser.add(affinefile);

    Option<string> warpfile("--warp", "", "Relative filename of warp (NB: Specify inverse warp, i.e. standard-to-native)", false, requires_argument);
    parser.add(warpfile);

    Option<bool> verbhelp("-h,--help", false, "Show help for this mode", false, no_argument);
    parser.add(verbhelp);

    int nextarg;
    try
    {
        nextarg = parser.parse_command_line(argc, argv, 0, true);
    }
    catch (X_OptionError &e)
    {
        exiterror(parser, e.what());
    }

    if (verbhelp.set())
    {
        parser.usage();
        return;
    }

    if (!parser.check_compulsory_arguments(true))
        exiterror(parser, "Not all compulsory arguments were set");

    if (nextarg == argc)
        exiterror(parser, "No subjects specified");

    if (maskfiles.set() && (affinefile.set() == warpfile.set()))
        exiterror(parser, "Either --affine or --warp should be set if standard space masks are used");

    std::vector<volume<float> > masks;
    if (maskfiles.set())
    {
        for (const auto &maskfile : splitarg(maskfiles.value()))
        {
            volume<float> mask;
            read_volume(mask, maskfile);
            masks.push_back(mask);
        }
    }

    std::vector<std::string> splitimagefiles = splitarg(imagefiles.value());
    
    for (; nextarg != argc; nextarg++)
    {
        std::string subj(argv[nextarg]);

        BOOST_LOG_TRIVIAL(info) << "Processing subject " << subj;
    
        std::vector<volume<float> > nativemasks;
        if (nativemaskfiles.set())
        {
            for (const auto &nativemaskfile : splitarg(nativemaskfiles.value()))
            {
                volume<float> nativemask;
                read_volume(nativemask, subj + "/" + nativemaskfile);
                nativemask.binarise(0.5);
                nativemasks.push_back(nativemask);
            }
        }

        std::vector<volume<float> > images;
        for (auto imagefile = splitimagefiles.cbegin(); imagefile != splitimagefiles.cend(); imagefile++)
        {
            volume<float> image;
            read_volume(image, subj + "/" + *imagefile);
            
            if (!images.empty()
                    && (image.xsize() != images[0].xsize()
                    || image.ysize() != images[0].ysize()
                    || image.zsize() != images[0].zsize()))
            {
                exiterror(parser, "All images should have the same dimensions.");
            }
            
            images.push_back(image);
        }

        if (!masks.empty())
        {
            boost::shared_ptr<const Transformation> transformation;
            if (affinefile.set())
                transformation = boost::make_shared<const AffineTransformation>(subj + "/" + affinefile.value());
            else
                transformation = boost::make_shared<const NonlinearTransformation>(subj + "/" + warpfile.value(), true);

            for (const auto &mask : masks)
            {
                volume<float> transformedmask;
                
                transformedmask = transformation->TransformVolume(mask, images[0]);
                transformedmask.binarise(0.5);

                for (const auto &image : images)
                    std::cout << image.percentile(0.5, transformedmask) << " ";
            }
        }

        if (!nativemasks.empty())
        {
            for (const auto &nativemask : nativemasks)
            {
                for (const auto &image : images)
                {
                    if (nativemask.xsize() != image.xsize()
                        || nativemask.ysize() != image.ysize()
                        || nativemask.zsize() != image.zsize())
                    {
                        exiterror(parser, "Images and native masks should have the same dimensions.");
                    }

                    std::cout << image.percentile(0.5, nativemask) << " ";
                }
            }
        }

        std::cout << std::endl;
    }
}

