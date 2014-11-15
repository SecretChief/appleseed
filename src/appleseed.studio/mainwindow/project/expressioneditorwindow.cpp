
//
// This source file is part of appleseed.
// Visit http://appleseedhq.net/ for additional information and resources.
//
// This software is released under the MIT license.
//
// Copyright (c) 2014 Marius Avram, The appleseedhq Organization
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
#include "expressioneditorwindow.h"

// UI definition headers.
#include "ui_expressioneditorwindow.h"

// appleseed.studio headers.
#include "mainwindow/project/tools.h"

// appleseed.renderer headers.
#include "renderer/global/globallogger.h"
#include "renderer/modeling/material/disneymaterial.h"
#include "renderer/modeling/project/project.h"

// appleseed.foundation headers.
#include "foundation/utility/foreach.h"

// appleseed.shared headers.
#include "application/application.h"

// SeExpr headers.
#pragma warning (push)
#pragma warning (disable : 4267)    // conversion from 'size_t' to 'int', possible loss of data
#include "SeExpression.h"
#include "SeExprEditor/SeExprEditor.h"
#include "SeExprEditor/SeExprEdBrowser.h"
#include "SeExprEditor/SeExprEdControlCollection.h"
#pragma warning (pop)

// Qt headers.
#include <QCloseEvent>
#include <QDesktopServices>
#include <QFileDialog>
#include <QHBoxLayout>
#include <QLabel>
#include <QPushButton>
#include <QScrollArea>
#include <QShortcut>
#include <QString>
#include <QVBoxLayout>

// Boost headers.
#include "boost/filesystem/operations.hpp"
#include "boost/filesystem/path.hpp"

// Standard headers.
#include <algorithm>
#include <fstream>
#include <sstream>

using namespace appleseed::shared;
using namespace boost;
using namespace foundation;
using namespace renderer;
using namespace std;

namespace appleseed {
namespace studio {

ExpressionEditorWindow::ExpressionEditorWindow(
    const Project&  project,
    const QString&  widget_name,
    const string&   expression,
    QWidget*        parent)
  : QWidget(parent)
  , m_project(project)
  , m_widget_name(widget_name)
  , m_ui(new Ui::ExpressionEditorWindow())
  , m_show_examples(false)
{
    m_ui->setupUi(this);
    setWindowFlags(Qt::Tool);
    setAttribute(Qt::WA_DeleteOnClose);
    QHBoxLayout* root_layout = new QHBoxLayout(m_ui->scrollarea);
    QVBoxLayout* left_layout = new QVBoxLayout();
    QVBoxLayout* right_layout = new QVBoxLayout();
    root_layout->addLayout(left_layout);
    root_layout->addLayout(right_layout);

    // Expression controls.
    SeExprEdControlCollection *controls = new SeExprEdControlCollection();
    QScrollArea* controls_scrollarea = new QScrollArea(this);
    controls_scrollarea->setObjectName("expression_controls");
    controls_scrollarea->setMinimumHeight(200);
    controls_scrollarea->setWidgetResizable(true);
    controls_scrollarea->setWidget(controls);
    left_layout->addWidget(controls_scrollarea);

    // Clear, Save, Load, Example buttons.
    QHBoxLayout* file_buttonbox = new QHBoxLayout();
    QPushButton* clear_button = new QPushButton("Clear");
    clear_button->setToolTip("Ctrl+N");
    QShortcut* clear_shortcut = new QShortcut(QKeySequence("Ctrl+N"), this);
    connect(clear_button, SIGNAL(clicked()), SLOT(slot_clear_expression()));
    connect(clear_shortcut, SIGNAL(activated()), SLOT(slot_clear_expression()));
    QPushButton* save_button = new QPushButton("Save");
    save_button->setToolTip("Ctrl+S");
    QShortcut* save_shortcut = new QShortcut(QKeySequence("Ctrl+S"), this);
    connect(save_button, SIGNAL(clicked()), SLOT(slot_save_script()));
    connect(save_shortcut, SIGNAL(activated()), SLOT(slot_save_script()));
    QPushButton* load_button = new QPushButton("Load");
    load_button->setToolTip("Ctrl+O");
    QShortcut* load_shortcut = new QShortcut(QKeySequence("Ctrl+O"), this);
    connect(load_button, SIGNAL(clicked()), SLOT(slot_load_script()));
    connect(load_shortcut, SIGNAL(activated()), SLOT(slot_load_script()));
    QPushButton *help_button = new QPushButton("Help");
    connect(help_button, SIGNAL(clicked()), SLOT(slot_show_help()));
    QPushButton* examples_button = new QPushButton("Examples");
    connect(examples_button, SIGNAL(clicked()), SLOT(slot_show_examples()));

    file_buttonbox->addWidget(clear_button);
    file_buttonbox->addWidget(save_button);
    file_buttonbox->addWidget(load_button);
    file_buttonbox->addWidget(help_button);
    file_buttonbox->addWidget(examples_button);
    left_layout->addLayout(file_buttonbox);

    m_editor = new SeExprEditor(this, controls);
    QTextEdit* text_edit = m_editor->findChild<QTextEdit*>("");
    text_edit->setObjectName("expression_editor");
    m_editor->setExpr(expression, true);
    left_layout->addWidget(m_editor);

    m_error = new QLabel("SeExpression has errors. View log for details.");
    m_error->setObjectName("error");
    m_error->hide();
    left_layout->addWidget(m_error);

    // Expression browser.
    m_browser = new SeExprEdBrowser(0, m_editor);
    const filesystem::path root_path(Application::get_root_path());
    const string scripts_path = (root_path / "seexpr").string();
    m_browser->addPath("Examples", scripts_path);
    m_browser->update();
    m_browser->hide();
    right_layout->addWidget(m_browser);

    m_ui->buttonbox_layout->addStretch(1);
    m_ui->buttonbox_layout->setStretch(0, 1);
    m_ui->buttonbox_layout->setStretch(1, 0);

    QPushButton* apply_button = m_ui->buttonbox->button(QDialogButtonBox::Apply);

    // Create connections.
    connect(m_ui->buttonbox, SIGNAL(accepted()), SLOT(slot_accept()));
    connect(apply_button, SIGNAL(clicked()), SLOT(slot_apply()));
    connect(m_ui->buttonbox, SIGNAL(rejected()), SLOT(slot_cancel()));
}

void ExpressionEditorWindow::apply_expression()
{
    const string expression = m_editor->getExpr();
    const DisneyParamExpression se_expression(expression.c_str());

    if (se_expression.is_valid())
    {
        m_error->hide();
        RENDERER_LOG_INFO("SeExpression successfully applied.");
        const QString q_expression = QString::fromStdString(expression);
        emit signal_expression_applied(m_widget_name, q_expression);
    }
    else
    {
        m_error->show();
        se_expression.report_error("SeExpression has errors");
    }
}

void ExpressionEditorWindow::slot_accept()
{
    apply_expression();
    close();
}

void ExpressionEditorWindow::slot_apply()
{
    apply_expression();
}

void ExpressionEditorWindow::slot_cancel()
{
    close();
}

void ExpressionEditorWindow::slot_clear_expression()
{
    m_editor->setExpr("");
}

void ExpressionEditorWindow::slot_save_script()
{
    QFileDialog::Options options;
    QString selected_filter;

    if (m_script_filepath.empty())
    {
        QString filepath =
            QFileDialog::getSaveFileName(
                this,
                "Save As...",
                QString::fromStdString(get_project_path()),
                "Expression Scripts (*.se)",
                &selected_filter,
                options);

        if (!filepath.isEmpty())
        {
            filepath = QDir::toNativeSeparators(filepath);
            m_script_filepath = filepath.toStdString();
        }
    }

    if (!m_script_filepath.empty())
    {
        ofstream script_file(m_script_filepath.c_str());
        if (!script_file.is_open())
        {
            show_warning_message_box(
                "Saving Error",
                "Failed to save the expression script file " + m_script_filepath + ".");
            return;
        }
        script_file << m_editor->getExpr();
        script_file.close();
    }
}

void ExpressionEditorWindow::slot_load_script()
{
    QFileDialog::Options options;
    QString selected_filter;

    QString filepath =
        QFileDialog::getOpenFileName(
            this,
            "Open...",
            QString::fromStdString(get_project_path()),
            "Expression Scripts (*.se);;All Files (*.*)",
            &selected_filter,
            options);

    if (!filepath.isEmpty())
    {
        filepath = QDir::toNativeSeparators(filepath);

        // Read script and set it as an expression.
        ifstream script_file(filepath.toStdString().c_str());
        if (!script_file.is_open())
        {
            show_warning_message_box(
                "Loading Error",
                "Failed to load the expression script file " + filepath.toStdString() + ".");
            return;
        }
        stringstream script_buffer;
        script_buffer << script_file.rdbuf();
        script_file.close();

        m_editor->setExpr(script_buffer.str());
        apply_expression();
    }
}

void ExpressionEditorWindow::slot_show_examples()
{
    m_show_examples = !m_show_examples;
    if (m_show_examples)
    {
        if (width() < 800)
            resize(800, height());
        m_ui->buttonbox_layout->setStretch(0, 1);
        m_ui->buttonbox_layout->setStretch(1, 1);
        m_browser->show();
    }
    else
    {
        if (width() > 400)
            resize(max(400, width() - 400), height());
        m_ui->buttonbox_layout->setStretch(1, 0);
        m_browser->hide();
    }
}

void ExpressionEditorWindow::slot_show_help()
{
    filesystem::path docs_path = Application::get_root_path();
    docs_path /= "documentation/seexpr/userdoc.html";
    QString docs_file = QString::fromStdString(docs_path.string());
    QDesktopServices::openUrl(QUrl::fromLocalFile(docs_file));
}

string ExpressionEditorWindow::get_project_path() const
{
    const filesystem::path project_root_path = filesystem::path(m_project.get_path()).parent_path();
    const filesystem::path file_path = absolute("script.se", project_root_path);
    const filesystem::path file_root_path = file_path.parent_path();

    return file_root_path.string();
}

void ExpressionEditorWindow::closeEvent(QCloseEvent *e)
{
    emit editor_closed();
    e->accept();
}

}       // namespace studio
}       // namespace appleseed

#include "mainwindow/project/moc_cpp_expressioneditorwindow.cxx"
