#!/usr/bin/env python
# coding: utf-8

# In[33]:


def equal_subsample(df, stimu, gt):
    temp = pd.concat([df[(df["genotype"]=="WT")&(df["stim"] == stimu)].sample(1250,random_state = seed),
                     df[(df["genotype"]=="KO")&(df["stim"] == stimu)].sample(1250,random_state = seed)], axis=0)
    temp = temp[temp["genotype"].isin(gt)]
    
    return temp


# In[41]:


### import libraries
import pandas as pd
from bokeh.models import ColorBar,HoverTool, ColumnDataSource
from bokeh.models.mappers import LinearColorMapper
from bokeh.palettes import Viridis
from bokeh.io import show, output_notebook, reset_output
from bokeh.plotting import Figure
from bokeh.io import curdoc
from bokeh.layouts import row, column
from bokeh.models.widgets import CheckboxGroup, RadioButtonGroup, Paragraph
seed = 42
umap_df = pd.concat([pd.read_csv("data/umap.csv",index_col=0),
                     pd.read_csv("data/umap_unstim.csv",index_col=0)],axis=0).reset_index(drop=True)


plot = Figure(title="Umap", x_axis_label ="umap_dim_1", y_axis_label ="umap_dim_2",
                plot_width=800, plot_height=800,toolbar_sticky=False, toolbar_location="below")

tooltips = [("Population","@label")]
plot.add_tools(HoverTool(tooltips=tooltips))

geno_title = Paragraph(text="Genotype: ", width=20, height=10)
genotypes = list(umap_df.genotype.unique())
types_checkbox = CheckboxGroup(labels=genotypes, active = [0,1])

ag_title = Paragraph(text="Antigen: ", width=20, height=10)
antigens = list(umap_df.columns[13:19])
ag_button = RadioButtonGroup(labels=antigens, active = 1)
active = ag_button.labels[ag_button.active]

stim_title = Paragraph(text="Treatment: ", width=20, height=10)
stim = ["stimulated", "Unstimulated"]
stim_button = RadioButtonGroup(labels=stim, active = 0)
stim_status = stim_button.labels[stim_button.active]

df = equal_subsample(umap_df, stim_status, [types_checkbox.labels[i] for i in types_checkbox.active])

source = ColumnDataSource({elm:df[elm].values for elm in df.columns})

mapper = LinearColorMapper(palette=Viridis[256], low=umap_df[active].min(), high=umap_df[active].max()/3)
cb = ColorBar(color_mapper=mapper, label_standoff=7,major_tick_line_color="black",margin=-10, title="MFI "+active,
             title_text_font_size = "8pt",title_text_font_style="bold",major_label_text_font_style="bold")
plot.add_layout(cb, 'right')
c = {'field': active, 'transform': mapper}
circ = plot.circle("umap_dim_1", "umap_dim_2", size=8, 
            source= source, fill_color=c,line_color="Black", line_width=0.5)

plot.xaxis.major_tick_line_color = None  # turn off x-axis major ticks
plot.xaxis.minor_tick_line_color = None  # turn off x-axis minor ticks

plot.yaxis.major_tick_line_color = None  # turn off y-axis major ticks
plot.yaxis.minor_tick_line_color = None  # turn off y-axis minor ticks

plot.xaxis.major_label_text_font_size = '0pt'  # turn off x-axis tick labels
plot.yaxis.major_label_text_font_size = '0pt'  # turn off y-axis tick labels

plot.xaxis.axis_label_text_font_size = '0pt'
plot.yaxis.axis_label_text_font_size = '0pt'

def update_plot():
    
    selected_types = [types_checkbox.labels[i] for i in types_checkbox.active]
    new_stim_status =stim_button.labels[stim_button.active]
    
    new_df = equal_subsample(umap_df, new_stim_status, selected_types)

    source.data = {elm:new_df[elm].values for elm in new_df.columns}
    new_active = ag_button.labels[ag_button.active]
    mapper.low = umap_df[new_active].min()
    mapper.high = umap_df[new_active].max()/3
    
    cb.title = "MFI "+new_active
    
    c_new = {'field': new_active, 'transform': mapper}
    circ.glyph.fill_color = c_new
    
active_changes = [types_checkbox, ag_button, stim_button]
for change in active_changes:
    change.on_change('active', lambda attr, old, new: update_plot())

layout = column(row(column(geno_title,types_checkbox), column(stim_title,stim_button)),row(column(ag_title,ag_button)),
                row(plot))
curdoc().add_root(layout)

