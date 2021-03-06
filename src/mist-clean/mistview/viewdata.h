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

#ifndef VIEWDATA_H
#define VIEWDATA_H

#include <vtkSmartPointer.h>
#include <vtkImageData.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <newimage.h>
#include <string>
#include <boost/shared_ptr.hpp>

class vtkLookupTable;

class ViewData
{
public:
    template <class T> ViewData(boost::shared_ptr<NEWIMAGE::volume<T> > vol, const std::string& displayname) :
        m_displayName(displayname)
    {
        m_imageData = NewimageToVtk(*vol);
    }

    vtkSmartPointer<vtkImageData> GetImageData();
    std::string GetDisplayName() const;
    double GetDisplayMin() const;
    double GetDisplayMax() const;
    void SetDisplayMin(const double& min);
    void SetDisplayMax(const double& max);
    vtkSmartPointer<vtkLookupTable> GetLut();
    void SetBinaryLut(double hue, double saturation);

private:
    vtkSmartPointer<vtkImageData> m_imageData;
    std::string m_displayName;
    vtkSmartPointer<vtkLookupTable> m_lut;

    static vtkSmartPointer<vtkLookupTable> BuildLut(double min, double max);
    static vtkSmartPointer<vtkLookupTable> BuildBinaryLut(double hue, double saturation);

    template <class T> vtkSmartPointer<vtkImageData> NewimageToVtk(const NEWIMAGE::volume<T>& vol)
    {
        vtkSmartPointer<vtkFloatArray> values = vtkSmartPointer<vtkFloatArray>::New();

        values->SetNumberOfComponents(1);
        for (int k = vol.minz(); k <= vol.maxz(); k++)
            for (int j = vol.miny(); j <= vol.maxy(); j++)
                for (int i = vol.minx(); i <= vol.maxx(); i++)
                   values->InsertNextTuple1(static_cast<float>(vol(i, j, k)));

        vtkSmartPointer<vtkImageData> imagedata = vtkSmartPointer<vtkImageData>::New();
        imagedata->SetExtent(0, vol.xsize() - 1, 0, vol.ysize() - 1, 0, vol.zsize() - 1);
        imagedata->SetSpacing(vol.xdim(), vol.ydim(), vol.zdim());
        imagedata->GetPointData()->SetScalars(values);

        m_lut = BuildLut(vol.robustmin(), vol.robustmax());

        return imagedata;
    }
};

#endif // VIEWDATA_H
