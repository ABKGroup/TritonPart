///////////////////////////////////////////////////////////////////////////////
// BSD 3-Clause License
//
// Copyright (c) 2019, The Regents of the University of California
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice, this
//   list of conditions and the following disclaimer.
//
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
//
// * Neither the name of the copyright holder nor the names of its
//   contributors may be used to endorse or promote products derived from
//   this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#pragma once

#include <tcl.h>

#include <QDockWidget>
#include <QPlainTextEdit>
#include <QPushButton>
#include <QSettings>
#include <QStringList>

#include "tclCmdInputWidget.h"
#include "utl/Logger.h"

namespace odb {
class dbDatabase;
}  // namespace odb

namespace gui {

// This shows a line edit to enter tcl commands and a
// text area that is used to show the commands and their
// results.  A command history is maintained and stored in
// settings across runs.  Up/down arrows scroll through the
// history as usual.  Qt itself provides editing features
// within the line edit.
class ScriptWidget : public QDockWidget
{
  Q_OBJECT

 public:
  ScriptWidget(QWidget* parent = nullptr);
  ~ScriptWidget();

  void readSettings(QSettings* settings);
  void writeSettings(QSettings* settings);

  void setLogger(utl::Logger* logger);

  void setupTcl(Tcl_Interp* interp,
                bool interactive,
                bool do_init_openroad,
                const std::function<void(void)>& post_or_init);

  void setFont(const QFont& font);

  void bufferOutputs(bool state);

 signals:
  // Commands might have effects that others need to know
  // (eg change placement of an instance requires a redraw)
  void commandExecuted(int return_code);
  void commandAboutToExecute();
  void executionPaused();

  // tcl exit has been initiated, want the gui to handle
  // shutdown
  void tclExiting();

  void addToOutput(const QString& text, const QColor& color);

 public slots:
  // Triggered when the user hits return in the line edit
  void executeCommand(const QString& command, bool echo = true);

  // Use to execute a command silently, ie. without echo or return.
  void executeSilentCommand(const QString& command);

 private slots:
  void outputChanged();

  void pause(int timeout);
  void unpause();

  void pauserClicked();

  void goBackHistory();
  void goForwardHistory();

  void updatePauseTimeout();

  void addTextToOutput(const QString& text, const QColor& color);

 protected:
  // required to ensure input command space it set to correct height
  void resizeEvent(QResizeEvent* event) override;

 private:
  int executeTclCommand(const QString& command);

  void triggerPauseCountDown(int timeout);

  void addCommandToOutput(const QString& cmd);
  void addTclResultToOutput(int return_code);
  void addReportToOutput(const QString& text);
  void addLogToOutput(const QString& text, const QColor& color);

  static int tclExitHandler(ClientData instance_data,
                            Tcl_Interp* interp,
                            int argc,
                            const char** argv);

  QPlainTextEdit* output_;
  TclCmdInputWidget* input_;
  QPushButton* pauser_;
  std::unique_ptr<QTimer> pause_timer_;
  Tcl_Interp* interp_;
  QStringList history_;
  QString history_buffer_last_;
  int historyPosition_;
  bool paused_;
  utl::Logger* logger_;

  bool buffer_outputs_;
  bool is_interactive_;

  // Logger sink
  template <typename Mutex>
  class GuiSink;
  std::shared_ptr<spdlog::sinks::sink> sink_;

  // maximum number of character to display in a log line
  const int max_output_line_length_ = 1000;

  const QColor cmd_msg_ = Qt::black;
  const QColor tcl_error_msg_ = Qt::red;
  const QColor tcl_ok_msg_ = Qt::blue;
  const QColor buffer_msg_ = QColor(0x30, 0x30, 0x30);
};

}  // namespace gui