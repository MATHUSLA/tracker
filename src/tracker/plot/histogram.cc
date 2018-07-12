/*
 * src/tracker/plot/histogram.cc
 *
 * Copyright 2018 Brandon Gomes
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <tracker/plot.hh>

#include <unordered_map>

#include <ROOT/TCanvas.h>
#include <ROOT/TH1D.h>
#include <ROOT/TFile.h>
#include <ROOT/TAxis.h>

namespace MATHUSLA { namespace TRACKER {

namespace plot { ///////////////////////////////////////////////////////////////////////////////

namespace { ////////////////////////////////////////////////////////////////////////////////////

//__Build TCanvas_______________________________________________________________________________
TCanvas* _build_TCanvas(const std::string& name,
                        const std::string& title) {
  return new TCanvas(name.c_str(), title.c_str(), 800, 800);
}
//----------------------------------------------------------------------------------------------

//__Build TH1D__________________________________________________________________________________
TH1D* _build_TH1D(const std::string& name,
                  const std::string& title,
                  const size_t bins,
                  const real min,
                  const real max) {
  return new TH1D(name.c_str(), title.c_str(), bins, min, max);
}
//----------------------------------------------------------------------------------------------

} /* anonymous namespace */ ////////////////////////////////////////////////////////////////////

//__Histogram Implementation Definition_________________________________________________________
struct histogram::impl {
  TCanvas* _canvas;
  TH1D* _hist;
  bool _has_updated = false;

  TAxis* x_axis() { return _hist->GetXaxis(); }
  const TAxis* x_axis() const { return _hist->GetXaxis(); }
  TAxis* y_axis() { return _hist->GetYaxis(); }
  const TAxis* y_axis() const { return _hist->GetYaxis(); }

  impl() = default;

  impl(const std::string& name,
       const std::string& title,
       const std::string& x_title,
       const std::string& y_title,
       const size_t bins,
       const real min,
       const real max)
    : _canvas(_build_TCanvas(name, title)),
      _hist(_build_TH1D(name, title, bins, min, max)) {
    _hist->SetDirectory(nullptr);
    x_axis()->SetTitle(x_title.c_str());
    y_axis()->SetTitle(y_title.c_str());
  }

  impl(const impl& other)
      : _canvas(_build_TCanvas(other._hist->GetName(), other._hist->GetTitle())),
        _hist(dynamic_cast<TH1D*>(other._hist->Clone())) {}

  impl(impl&& other) noexcept = default;

  impl& operator=(const impl& other) {
    if (this != &other) {
      _canvas = dynamic_cast<TCanvas*>(other._canvas->Clone());
      _hist = dynamic_cast<TH1D*>(other._hist->Clone());
    }
    return *this;
  }

  impl& operator=(impl&& other) noexcept = default;
  ~impl() = default;
};
//----------------------------------------------------------------------------------------------

//__Histogram Empty Constructor_________________________________________________________________
histogram::histogram() : _impl(std::make_unique<impl>()) {}
//----------------------------------------------------------------------------------------------

//__Histogram Partial Constructor_______________________________________________________________
histogram::histogram(const histogram::name_type& name,
                     const size_t bins,
                     const real min,
                     const real max)
    : histogram(name, name, "X", "Y", bins, min, max) {}
//----------------------------------------------------------------------------------------------

//__Histogram Partial Constructor_______________________________________________________________
histogram::histogram(const histogram::name_type& name,
                     const std::string& title,
                     const size_t bins,
                     const real min,
                     const real max)
    : histogram(name, title, "X", "Y", bins, min, max) {}
//----------------------------------------------------------------------------------------------

//__Histogram Full Constructor__________________________________________________________________
histogram::histogram(const histogram::name_type& name,
                     const std::string& title,
                     const std::string& x_title,
                     const std::string& y_title,
                     const size_t bins,
                     const real min,
                     const real max)
    : _impl(std::make_unique<impl>(name, title, x_title, y_title, bins, min, max)) {}
//----------------------------------------------------------------------------------------------

//__Histogram Copy Constructor__________________________________________________________________
histogram::histogram(const histogram& other) : _impl(new impl(*other._impl)) {}
//----------------------------------------------------------------------------------------------

//__Histogram Move Constructor__________________________________________________________________
histogram::histogram(histogram&& other) noexcept = default;
//----------------------------------------------------------------------------------------------

//__Histogram Copy Assignment___________________________________________________________________
histogram& histogram::operator=(const histogram& other) {
  if (this != &other)
    _impl.reset(new impl(*other._impl));
  return *this;
}
//----------------------------------------------------------------------------------------------

//__Histogram Move Assignment___________________________________________________________________
histogram& histogram::operator=(histogram&& other) noexcept = default;
//----------------------------------------------------------------------------------------------

//__Histogram Destructor________________________________________________________________________
histogram::~histogram() = default;
//----------------------------------------------------------------------------------------------

//__Get Histogram Name__________________________________________________________________________
const histogram::name_type histogram::name() const {
  return _impl->_hist->GetName();
}
//----------------------------------------------------------------------------------------------

//__Get Histogram Title_________________________________________________________________________
const std::string histogram::title() const {
  return _impl->_hist->GetTitle();
}
//----------------------------------------------------------------------------------------------

//__Get Histogram X-Axis Title__________________________________________________________________
const std::string histogram::x_title() const {
  return _impl->x_axis()->GetTitle();
}
//----------------------------------------------------------------------------------------------

//__Get Histogram Y-Axis Title__________________________________________________________________
const std::string histogram::y_title() const {
  return _impl->y_axis()->GetTitle();
}
//----------------------------------------------------------------------------------------------

//__Set Histogram Name__________________________________________________________________________
void histogram::name(const histogram::name_type& name) {
  _impl->_hist->SetName(name.c_str());
}
//----------------------------------------------------------------------------------------------

//__Set Histogram Title_________________________________________________________________________
void histogram::title(const std::string& title) {
  _impl->_hist->SetTitle(title.c_str());
}
//----------------------------------------------------------------------------------------------

//__Set Histogram X-Axis Title__________________________________________________________________
void histogram::x_title(const std::string& x_title) {
  _impl->x_axis()->SetTitle(x_title.c_str());
}
//----------------------------------------------------------------------------------------------

//__Set Histogram Y-Axis Title__________________________________________________________________
void histogram::y_title(const std::string& y_title) {
  _impl->y_axis()->SetTitle(y_title.c_str());
}
//----------------------------------------------------------------------------------------------

//__Check If Histogram is Empty_________________________________________________________________
bool histogram::empty() const {
  return size() == 0;
}
//----------------------------------------------------------------------------------------------

//__Size of Histogram___________________________________________________________________________
size_t histogram::size() const {
  return _impl->_hist->GetEntries();
}
//----------------------------------------------------------------------------------------------

//__Least Bin of Histogram______________________________________________________________________
real histogram::min_x() const {
  const auto& hist = _impl->_hist;
  return hist->GetBinLowEdge(hist->GetMinimumBin());
}
//----------------------------------------------------------------------------------------------

//__Greatest Bin of Histogram___________________________________________________________________
real histogram::max_x() const {
  const auto& hist = _impl->_hist;
  const auto max_bin = hist->GetMaximumBin();
  return hist->GetBinLowEdge(max_bin) + hist->GetBinWidth(max_bin);
}
//----------------------------------------------------------------------------------------------

//__Average of Histogram________________________________________________________________________
real histogram::mean() const {
  return _impl->_hist->GetMean();
}
//----------------------------------------------------------------------------------------------

//__Value at Given Bin__________________________________________________________________________
real histogram::bin_value(const size_t index) const {
  return _impl->_hist->GetBinContent(index);
}
//----------------------------------------------------------------------------------------------

//__Value at Given Point________________________________________________________________________
real histogram::value(const real point) const {
  return bin_value(_impl->_hist->FindBin(point));
}
//----------------------------------------------------------------------------------------------

//__Insert Point into Histogram_________________________________________________________________
size_t histogram::insert(const real point) {
  return _impl->_hist->Fill(point);
}
//----------------------------------------------------------------------------------------------

//__Increment Bin_______________________________________________________________________________
void histogram::increment(const size_t index,
                          const real weight) {
  _impl->_hist->AddBinContent(index, weight);
}
//----------------------------------------------------------------------------------------------

//__Scale Histogram_____________________________________________________________________________
void histogram::scale(const real weight) {
  _impl->_hist->Scale(weight);
}
//----------------------------------------------------------------------------------------------

//__Draw Histogram to Canvas____________________________________________________________________
void histogram::draw() {
  // TODO: think about plot::is_on uses
  if (plot::is_on()) {
    _impl->_canvas->cd();
    _impl->_hist->Draw("HIST");
    _impl->_canvas->Modified();
    _impl->_canvas->Update();
    _impl->_has_updated = true;
  }
}
//----------------------------------------------------------------------------------------------

//__Clear Histogram and Canvas__________________________________________________________________
void histogram::clear() {
  if (plot::is_on() && _impl->_has_updated) {
    _impl->_canvas->cd();
    _impl->_canvas->Clear();
    _impl->_canvas->Modified();
    _impl->_canvas->Update();
    _impl->_hist->Reset("ICES");
    _impl->_has_updated = false;
  }
}
//----------------------------------------------------------------------------------------------

//__Save Histogram to File______________________________________________________________________
bool histogram::save(const std::string& path) const {
  TFile file(path.c_str(), "UPDATE");
  if (!file.IsZombie()) {
    file.cd();
    file.WriteTObject(_impl->_hist->Clone());
    file.Close();
    return true;
  }
  return false;
}
//----------------------------------------------------------------------------------------------

//__Load Histogram From File____________________________________________________________________
histogram histogram::load(const std::string& path,
                          const histogram::name_type& name) {
  histogram out;
  TFile file(path.c_str(), "READ");
  if (!file.IsZombie()) {
    TH1D* test = nullptr;
    file.GetObject(name.c_str(), test);
    if (test) {
      out._impl->_hist = test;
      out._impl->_canvas = _build_TCanvas(test->GetName(), test->GetTitle());
    }
    file.Close();
  }
  return out;
}
//----------------------------------------------------------------------------------------------

} /* namespace plot */  ////////////////////////////////////////////////////////////////////////

} } /* namespace MATHUSLA::TRACKER */