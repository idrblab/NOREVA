/**
 * Copyright 1997-2021 <a href="mailto:zhanghn@zju.edu.cn">Zhang Hongning</a>.
 * 
 * Modified at 2020-02-04
 * Licensed under the Apache License, Version 2.0 (thie "License");
 * You may not use this file except in compliance with the license.
 * You may obtain a copy of the License at
 * 
 *       http://wwww.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language govering permissions and
 * limitations under the License.
 */
package biojar.function.graphics;

import static biojar.application.SettingFrame.getDefaultDelimiter;
import static biojar.function.graphics.ImageHelp.getSystemFont;

import java.awt.Color;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics2D;
import java.io.File;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.LineNumberReader;

import java.util.ArrayList;

import biojar.function.GeneralMethod;

/**
 * 绘制环状柱形图的业务实现类
 * @version 2.0
 * @since 14 2021-02-04
 * @author <a href="mailto:zhanghn@zju.edu.cn">Zhang Hongning</a>
 */
public class CircularHistogram {
	/**
	 * 背景颜色，默认为白色（255,255,255）
	 */
	private Color bgColor = Color.WHITE;
	/**
	 * 字体颜色，默认为黑色（0,0,0）
	 */
	private Color fontColor = Color.BLACK;
	/**
	 * 总角度为340°, 但应满足totalAngle + 2*angle 不超过360，否则会对totalAngle自动调整
	 */
	private double totalAngle = 340;
	/**
	 * 颜色对象数组，用于对不同标准的上色
	 */
	private Color[] barColorSet = {
		new Color(128,0,128),//紫色 #800080，最外圈颜色
		new Color(251,188,5),//黄色 #FBBC05，次外圈颜色
		new Color(66,133,244),//蓝色 #4285F4
		new Color(234,67,53),//红色 #EA4335
		new Color(189,183,107),//深卡其布 #BDB76B
		new Color(107,142,35),//深绿 #6B8E23
		new Color(135,206,250),//淡蓝色 #87CEFA
		new Color(186,186,186)//灰色 #BABABA
	};
	/**
	 * 图例文本字符串组，第一个元素为图例标题，其余为各项文本
	 */
	private String[] legendTextSet = null;
	/**
	 * 柱形最大长度所表示的特征值
	 */
	private double maxValue = 40;
	/**
	 * 2D绘图对象，用于绘制图形
	 */
	private Graphics2D graphics = null;
	/**
	 * 输入数据存储ArrayList，元素类型为Object[]
	 * Object[]对象length固定为6，依次为String, int, int, int, int, int
	 */
	private ArrayList<Object[]> criteriaList = new ArrayList<>();
	/**
	 * 字体指标，用于对字体绘图进行定位
	 */
	private FontMetrics metrics;
	/**
	 * 输出文件名
	 */
	private String outfilename = null;
	/**
	 * 对jpg图像进行dpi设置，单位为像素/英寸
	 */
	@SuppressWarnings("unused")
	private int dpi = 360;
	/**
	 * 返回输出文件对象
	 * @return 输出文件对象
	 */
	public File getOutputFile() {
		if (outfilename != null) return new File(outfilename);
		return null;
	}
	/**
	 * 设置jpg图像的输出dpi，单位为像素/英寸，包括纵向和横向dpi
	 * @param dp dpi值，要求为正整数
	 */
	public void setJpegDpi(Integer dp) {
		if (dp <=0) {
			System.err.println("dpi should be more than 0.");
			return;
		}
		dpi = dp;
	}
	/**
	 * 设置图例文本
	 * @param args 新的图例文本字符串组
	 */
	public void setLegendTextSet(String[] args) {
		if (legendTextSet == null) legendTextSet = args;
	}
	/**
	 * 设置图例文件
	 * @param length 输入未定义图例时自定设置图例文本
	 */
	public void setLegendTextSet (Integer length) {
		if (legendTextSet == null) {
			legendTextSet = new String[length];
			legendTextSet[0] = "Label";
			for (int i = 1; i < length; i++) {
				legendTextSet[i] = "Label " + i;
			}
		}
	}
	/**
	 * 设置柱形最大长度所表示的特征值，至少是特征值的最大值，当设定值小于输入的最大值，会自动变为最大值而使得设置失效
	 * @param max 柱形最大长度所表示的新特征值，要求是正数
	 */
	public void setMaxValue(Double max) {
		if (max <=0) {
			System.err.println("Max value should be more than 0.0.");
			return;
		}
		maxValue = max;
	}
	/**
	 * 自定义设置总角度，单位为度
	 * @param ag 新的总角度值，要求在0.0-360.0之间，但应满足totalAngle + 2*angle 不超过360，
	 * 否则会对totalAngle自动调整，其中angle指的是旋转步长。这在cutoff值较小或输入数量较小时会发生
	 */
	public void setTotalAngle(Double ag) {
		if (ag<=0||ag>360) {
			System.err.println("Total angle should between 0.0 (not included) and 360.0 degree.");
			return;
		}
		totalAngle = ag;
	}
	/**
	 * 设置背景颜色
	 * @param color 十六进制颜色
	 */
	public void setBackgroundColor(String color) {
		bgColor = Color.decode(color);
	}
	/**
	 * 设置字体颜色
	 * @param color 十六进制颜色
	 */
	public void setFontColor(String color) {
		fontColor = Color.decode(color);
	}
	/**
	 * 对图配色进行自定义设置
	 * @param colorSet 输入的颜色字符串组，格式为#CC00FF
	 */
	public void setBarColorSet(String[] colorSet) {
		try {
			for (int i=0; i < barColorSet.length & i < colorSet.length; i++)
				barColorSet[i] = Color.decode(colorSet[i]);
		} catch (NumberFormatException ex) {
			System.err.println("Ilegal hexadecimal string! it should be like \"#CC00FF\"");
		}
	}
	/**
	 * 对图颜色进行自定义配置
	 * @param colorSet 输入的颜色组
	 */
	public void setBarColorSet(Color[] colorSet) {
		barColorSet = colorSet;
	}
	/**
	 * 绘制文本标签
	 * @param label 标签内容
	 * @param locat_x 标签定位横坐标（从左上角为0）
	 * @param locat_y 标签定位纵坐标（从左上角为0）
	 * @param h_mode 水平对齐模式符，m为水平居中对齐，l为水平左对齐，r为水平右对齐
	 * @param v_mode 垂直对齐模式符，m为垂直居中对齐，u为顶端对齐，d为底端对齐
	 * @param FontColor 文本标签字体颜色
	 * @throws Exception 异常
	 */
	private void drawSimpleLabel (String label, double rotateCenter_x, double rotateCenter_y, double rotateDegree, double rotateR, String h_mode, String v_mode, Color FontColor) throws Exception {
		float baseline_x;
		float baseline_y;
		int locat_x = (int) (rotateCenter_x - rotateR), locat_y = (int) rotateCenter_y;
		if (h_mode.equals("l")) {
			locat_x = (int) (rotateCenter_x + rotateR);
		}
		switch	 (v_mode) {
			case "m":{baseline_y =locat_y*1f - metrics.getHeight()/2f + metrics.getAscent();break;}
			case "u":{baseline_y =locat_y + metrics.getAscent();break;}
			case "d":{baseline_y =locat_y - metrics.getHeight() + metrics.getAscent();break;}
			default: throw new Exception("Ilegal mode symbol: "+v_mode);
		}
		switch (h_mode) {
			case "m":{baseline_x =locat_x*1f - metrics.stringWidth(label)/2f;break;}
			case "l":{baseline_x =locat_x;break;}
			case "r": {baseline_x =locat_x - metrics.stringWidth(label);break;}
			default: throw new Exception("Ilegal mode symbol: "+h_mode);
		}
		graphics.setColor(FontColor);//设置字体颜色
		graphics.rotate(rotateDegree*Math.PI/180, rotateCenter_x, rotateCenter_y);
		graphics.drawString(label, baseline_x, baseline_y);
		graphics.rotate(-rotateDegree*Math.PI/180, rotateCenter_x, rotateCenter_y);
	}
	/**
	 * 绘制简单旋转矩形
	 * @param rotateCenter_x 旋转中心横坐标
	 * @param rotateCenter_y 旋转中心纵坐标
	 * @param rotateDegree 旋转角度，单位为度，正数为顺时针旋转
	 * @param rotateR 旋转半径
	 * @param width 矩形宽度
	 * @param length 矩形长度
	 * @param color 矩形填充颜色
	 */
	private void drawSimpleBar(double rotateCenter_x, double rotateCenter_y, double rotateDegree, double rotateR, int width, int length, Color color) {
		graphics.setColor(color);
		if (width%2!=0) width++;
		int baseX = (int) (rotateCenter_x - width/2);
		int baseY = (int) (rotateCenter_y + rotateR);
		graphics.rotate(rotateDegree*Math.PI/180, rotateCenter_x, rotateCenter_y);//正数为顺时针转形状，也就是逆时针转画布, 旋转一次画一次
		graphics.fillRect(baseX, baseY, width, length);
		graphics.rotate(-rotateDegree*Math.PI/180, rotateCenter_x, rotateCenter_y);//因此转回来，方便统计总角度
	}
	/**
	 * 对数组进行排序，仅对内部数组使用
	 * @param arraylist 待排序数组
	 * @param order boolean值，true为升序，false为降序
	 * @return 排序后产生的ArrayList
	 */
	private static ArrayList<Object[]> sorted(ArrayList<Object[]> arraylist, boolean order) {
		ArrayList<Object[]> res = new ArrayList<Object[]>();
		@SuppressWarnings("unchecked")
		ArrayList<Object[]> list = (ArrayList<Object[]>) arraylist.clone();
		int length = list.size();
		if (length > 0) {
			for (int index = 0;index<length;index++) {
				Object[] value = list.get(0);
				for (Object[] e: list) {
					if (order) {
						if((double) value[1] >(double) e[1]) value = e;
					} else {
						if((double) value[1] <(double) e[1]) value = e;
					}
				}
				res.add(value);
				list.remove(list.indexOf(value));
			}
		}
		return res;
	}
	/**
	 * 加载数据
	 * @param inputfile 输入文件名
	 * @param hasTitle 是否包含标题
	 * @param cutoff 筛选值，取前cutoff名保留
	 * @throws FileNotFoundException 文件缺失异常
	 * @throws IOException 输入输出异常
	 */
	public void loadData(String inputfile, Boolean hasTitle, Integer cutoff) throws FileNotFoundException, IOException {
		if (cutoff <=1) {
			System.err.println("cutoff should be more than one.");
			cutoff = 100;
		}
		try (LineNumberReader lnr = GeneralMethod.BufferRead(inputfile)) {
			String line;
			while ((line = lnr.readLine()) != null) {
				if (hasTitle && lnr.getLineNumber() == 1) {
					CircularHistogram.this.setLegendTextSet(line.split(getDefaultDelimiter()));//将标题加载为图例文本
					continue;
				}
				String[] arr = line.split(getDefaultDelimiter());
				setLegendTextSet(arr.length);
				String label = arr[0];
				double Ct = 0.0;
				Object[] criteria = new Object[arr.length + 1];
				criteria[0] = label;
				for (int i =1; i < arr.length; i++) {
					double curr = Double.parseDouble(arr[i]);
					Ct += curr;
					criteria[i + 1] = curr;
				}
				criteria[1] = Ct;
				criteriaList.add(criteria);
			}
		}
		criteriaList = sorted(criteriaList, false);
		if (cutoff > criteriaList.size()) {
			System.err.println("cutoff value is more than all data and will be reset as "+ criteriaList.size());
			cutoff = criteriaList.size();
		}
		ArrayList <Object[]> selectList = new ArrayList<>();
		for (int index=0; index < cutoff && index<criteriaList.size(); index++) selectList.add(criteriaList.get(index));
		criteriaList = selectList;//销毁完整数组，只保留前几个
	}
	/**
	 * 从ArrayList java对象中装载数据，要求元素类型为字符串组
	 * @param list ArrayList对象，元素为字符串组，组内元素分别为label和属性值
	 * @param args 图例文本，第一个元素为图例标题，其余为图例各项文本
	 * @param cutoff 筛选值，取前cutoff名保留
	 */
	public void loadData(ArrayList<String[]> list, String[] args, Integer cutoff) {
		if (cutoff <=1) {
			System.err.println("cutoff should be more than one.");
			cutoff = 100;
		}
		CircularHistogram.this.setLegendTextSet(args);
		for (String[] arr: list) {
			String label = arr[0];
			double Ct = 0.0;
			Object[] criteria = new Object[arr.length + 1];
			criteria[0] = label;
			for (int i =1; i < arr.length; i++) {
				double curr = Double.parseDouble(arr[i]);
				Ct += curr;
				criteria[i + 1] = curr;
			}
			criteria[1] = Ct;
			criteriaList.add(criteria);
		}
		criteriaList = sorted(criteriaList, false);
		if (cutoff > criteriaList.size()) {
			System.err.println("cutoff value is more than all data and will be reset as "+ criteriaList.size());
			cutoff = criteriaList.size();
		}
		ArrayList <Object[]> selectList = new ArrayList<>();
		for (int index=0; index < cutoff; index++) selectList.add(criteriaList.get(index));
		criteriaList = selectList;
	}
	/**
	 * 对Graphics2D对象进行绘制环形柱状图，Grraphics2D类不同子类绘制结果类型不同
	 * @param r 绘制图最内环半径，为基础半径
	 * @param fontSize 标签字体字号
	 * @param width 画布宽度
	 * @param height 画布高度
	 * @param angle 旋转绘制角度步长，单位为度，正数为顺时针旋转
	 * @param type 输出文件类型，支持jpg/eps/pdf
	 * @throws Exception 相关异常
	 */
	private void drawFigure(int baseR, int fontSize, int barWidth, int width, int height, double angle) throws Exception {
		int center_x = width-height/2;
		int center_y = height/2;
		double baseAngle = 450-totalAngle;
		/*
		设置图片背景色为白色
		 */
		graphics.setBackground(bgColor);
		graphics.clearRect(0, 0, width, height);
		/*
		绘制四层统计图
		*/
		int r = baseR;
		Object[] obj = criteriaList.get(0);
		for (int i = 0; i < obj.length - 2; i++) {
			for (int index = 0; index < criteriaList.size(); index++) {
				obj = criteriaList.get(index);
				double currentValue = (double) obj[obj.length - 1 - i];
				int length = (int) (currentValue*baseR/maxValue);
				if (length > baseR) {
					System.err.println(length + "\t" + currentValue);
					length = baseR;
				}
				int colorindex =obj.length - 3 - i;//最大圈颜色用索引最小值
				while (colorindex >= barColorSet.length) colorindex -= barColorSet.length;
				drawSimpleBar(center_x,
					center_y,
					baseAngle+index*angle,
					r,
					barWidth,
					length,
					barColorSet[colorindex]
				);
			}
			r += baseR;
		}
		/*
		绘制中心圆圈，消除棱角
		*/
		int ir = (int) (Math.sqrt(1.0*baseR*baseR+barWidth*barWidth/4.0) + 0.5);
		graphics.setColor(bgColor);
		graphics.fillOval(center_x- ir, center_y - ir,2*(int)ir, 2*ir);
		/*
		绘制外围标签
		*/
		r += baseR/10;//定义外围标签与图间距为r/50
		for (int index = 0; index < criteriaList.size(); index++) {
			double currentAngle = 360-totalAngle + index*angle;
			drawSimpleLabel(
				(String)criteriaList.get(index)[0],
				center_x,
				center_y,
				currentAngle,//(currentAngle <= 90 || currentAngle>=270)?currentAngle:currentAngle - 180,
				r,
				"r",//(currentAngle <= 90 || currentAngle>=270)?"r":"l",
				"m",
				fontColor
			);
		}
		/*
		添加图注
		*/
		int legendX = width/25;
		int legendY = height/15;
		fontSize = height/50;
		Font legendFont = new Font(null, Font.BOLD, fontSize);//用默认字体
		graphics.setFont(legendFont);
		metrics = graphics.getFontMetrics();//字体换，标尺属性也得更新
		drawSimpleLabel(legendTextSet[0], legendX, legendY, 0, 0, "l", "m", fontColor);
		graphics.setFont(new Font(null, Font.BOLD, fontSize*2/3));
		metrics = graphics.getFontMetrics();
		for (int i=0; i < obj.length - 2; i++) {
			int colorIndex = i;
			while (colorIndex >= barColorSet.length) colorIndex -= barColorSet.length;
			drawSimpleBar(legendX, legendY+(i + 1) *fontSize*5/4, -90, 0, fontSize, fontSize, barColorSet[colorIndex]);
			drawSimpleLabel(legendTextSet[i + 1], legendX + fontSize*3/2, legendY+(i + 1) *fontSize*5/4, 0, 0, "l", "m", fontColor);
		}
	}
	/**
	 * 绘制环形柱状图
	 * @param outputfileformat 输出文件基础名，与维度，类型一起构成完整文件名
	 * @param autoSize 是否自动调整尺寸，仅对jpg格式有效
	 * @param type 输出格式，暂且支持jpg/eps/pdf
	 * @throws IOException 输入输出异常
	 * @throws Exception 其他异常
	 */
	public void drawFigure(String outputfileformat, Boolean autoSize, String type) throws IOException, Exception {
		int cutoff = criteriaList.size();
		double angle = totalAngle/(cutoff - 1);
		if (2*angle + totalAngle > 360)
			System.err.println("Current total angle will be adjusted automatically because of \"2×step angle + total angle > 360\".");
		while (2*angle + totalAngle > 360) {//对总角度进行自适应，使其能够直观看到起止位置
			totalAngle = totalAngle - 10;
			angle = totalAngle/(criteriaList.size() - 1);
		}
		int fontSize=56;
		int barWidth = fontSize/2;
		double widthPerHeight = 4.0/3;
		int heightPerR =  2*criteriaList.get(0).length +3;
		int r= (int) Math.round(barWidth*180/(1.5*angle*Math.PI));
		int height = heightPerR *r;
		int width = (int) (height*widthPerHeight);
		if (width > 23000 || !autoSize) {
			width = 23000;
			height = (int) (width/widthPerHeight);
			r = height/heightPerR;
			barWidth = (int) Math.round(Math.PI/180*angle*r*1.5);
			fontSize = barWidth*2;
		}
		type = type.toLowerCase();
		if (!type.equals("jpg")) {
			width = 2000;
			height = (int) (width/widthPerHeight);
			r = height/heightPerR;
			barWidth = (int) Math.round(Math.PI/180*angle*r*1.5);
			fontSize = barWidth*2;
		}
		if (barWidth %2 !=0) barWidth++;
		if (fontSize < 1) fontSize = 1;//保证最小字号
		if (fontSize > 20) fontSize = 20;//防止字号过大溢出
		if (barWidth < 1) barWidth = 1;//保证至少一像素宽度
		/*
		实例化当前类，配置字体
		*/
		String fontName = "Courier New";
		String[] fontlist = getSystemFont();
		boolean fontExist = false;
		for (String name: fontlist) if (fontName.equals(name)) fontExist = true;
		if (!fontExist) fontName = null;//字体不存在时用默认字体
		Font awtFont = new Font(fontName, Font.BOLD, fontSize);
		outfilename = String.format(outputfileformat, cutoff, type);
		CreateGraphics cg = new CreateGraphics(width, height, type, outfilename);
		graphics = cg.getGraphics2D();
		graphics.setFont(awtFont);
		metrics = graphics.getFontMetrics();
		drawFigure(r, fontSize, barWidth, width, height, angle);
		cg.saveToFlie();
	} 
}