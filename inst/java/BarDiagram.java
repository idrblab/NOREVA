package BarFigureDraw;

import BarFigureDraw.GeneralMethod;
import static BarFigureDraw.ImageHelp.saveAsJPEG;
import static BarFigureDraw.ImageHelp.getSystemFont;
import static BarFigureDraw.ImageHelp.saveAsEps;
import static BarFigureDraw.ImageHelp.getFontFileName;

import com.lowagie.text.Document;
import com.lowagie.text.Rectangle;
import com.lowagie.text.pdf.BaseFont;
import com.lowagie.text.pdf.FontMapper;
import com.lowagie.text.pdf.PdfContentByte;
import com.lowagie.text.pdf.PdfWriter;

import org.sourceforge.jlibeps.epsgraphics.EpsGraphics2D;

import java.awt.Color;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.image.BufferedImage;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.LineNumberReader;

import java.util.ArrayList;
import java.util.logging.Level;
import java.util.logging.Logger;
/**
.*.@file		BarDiagram.java
 * @author		Zhang Hongning
.* @email		zhang2016@zju.edu.cn
 * @version	2.0
 * @created	2020-8-28 11:09:21
.* @decription	about this java project
 */

public class BarDiagram {
	/**
	 * 总角度为默认340°
	 */
	private double totalAngle = 340;
	/**
	 * 颜色对象数组，用于对不同标准的上色
	 */
	private Color[] colorReflect = {
		new Color(234,67,53),//红色 #EA4335
		new Color(66,133,244),//蓝色 #4285F4
		new Color(251,188,5),//黄色 #FBBC05
		new Color(128,0,128)//紫色 #800080
	};
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
	 * 对jpg图像进行dpi设置，单位为像素/英寸
	 */
	private int dpi = 360;
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
	 * @param ag 新的总角度值，要求在0.0-360.0之间
	 */
	public void setTotalAngle(Double ag) {
		if (ag<=0||ag>360) {
			System.err.println("Total angle should between 0.0 (not included) and 360.0 degree.");
			return;
		}
		totalAngle = ag;
	}
	/**
	 * 对图配色进行自定义设置
	 * @param colorSet 输入的颜色字符串组，格式为#CC00FF
	 */
	public void setColorSet(String[] colorSet) {
		try {
			for (int i=0; i < colorReflect.length & i < colorSet.length; i++)
				colorReflect[i] = Color.decode(colorSet[i]);
		} catch (NumberFormatException ex) {
			System.err.println("Ilegal hexadecimal string! it should be like \"#CC00FF\"");
		}
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
		int baseline_x;
		int baseline_y;
		int locat_x = (int) (rotateCenter_x - rotateR), locat_y = (int) rotateCenter_y;
		if (h_mode.equals("l")) {
			locat_x = (int) (rotateCenter_x + rotateR);
		}
		switch	 (v_mode) {
			case "m":{baseline_y =locat_y - metrics.getHeight()/2 + metrics.getAscent();break;}
			case "u":{baseline_y =locat_y + metrics.getAscent();break;}
			case "d":{baseline_y =locat_y - metrics.getHeight() + metrics.getAscent();break;}
			default: throw new Exception("Ilegal mode symbol: "+v_mode);
		}
		switch (h_mode) {
			case "m":{baseline_x =locat_x - metrics.stringWidth(label)/2;break;}
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
	private static ArrayList sorted(ArrayList<Object[]> arraylist, boolean order) {
		ArrayList<Object[]> res = new ArrayList();
		ArrayList<Object[]> list = (ArrayList) arraylist.clone();
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
	private void loadData(String inputfile, boolean hasTitle, int cutoff) throws FileNotFoundException, IOException {
		try (LineNumberReader lnr = GeneralMethod.BufferRead(inputfile)) {
			String line;
			while ((line = lnr.readLine()) != null) {
				if (hasTitle && lnr.getLineNumber() == 1) continue;
				String[] arr = line.split("\t");
				String label = arr[0];
				double Ca = Double.parseDouble(arr[1]);
				double Cb = Double.parseDouble(arr[2]);
				double Cc = Double.parseDouble(arr[3]);
				double Cd = Double.parseDouble(arr[4]);
				double Ct = Ca + Cb + Cc + Cd;
				Object[] criteria = {label, Ct, Ca, Cb, Cc, Cd};
				criteriaList.add(criteria);
			}
		}
		criteriaList = sorted(criteriaList, false);
		ArrayList <Object[]> selectList = new ArrayList<>();
		for (int index=0; index < cutoff && index<criteriaList.size(); index++) selectList.add(criteriaList.get(index));
		criteriaList = selectList;//销毁完整数组，只保留前几个
	}
	/**
	 * 对Graphics2D对象进行绘制散射柱状图，Grraphics2D类不同子类绘制结果类型不同
	 * @param r 绘制图最内环半径，为基础半径
	 * @param fontSize 标签字体字号
	 * @param width 画布宽度
	 * @param height 画布高度
	 * @param angle 旋转绘制角度步长，单位为度，正数为顺时针旋转
	 * @param type 输出文件类型，支持jpg/eps/pdf
	 * @throws Exception 相关异常
	 */
	private void drawFigure(int baseR, int fontSize, int barWidth, int width, int height, double angle) throws Exception {
		Color BgColor = Color.WHITE;//背景色
		int center_x = width-height/2;
		int center_y = height/2;
		double baseAngle = 450-totalAngle;
		/*
		设置图片背景色为白色
		 */
		graphics.setBackground(BgColor);
		graphics.clearRect(0, 0, width, height);
		/*
		绘制四层统计图
		*/
		int r = baseR;
		for (int i = 0; i < 4; i++) {
			for (int index = 0; index < criteriaList.size(); index++) {
				double currentValue = (double) criteriaList.get(index)[5-i];
				if (currentValue > maxValue) maxValue = currentValue;
				int length = (int) (currentValue*baseR/maxValue);
				drawSimpleBar(
					center_x,
					center_y,
					baseAngle+index*angle,
					r,
					barWidth,
					length,
					colorReflect[i]
				);
			}
			r += baseR;
		}
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
				(currentAngle <= 90 || currentAngle>=270)?currentAngle:currentAngle - 180,
				r,
				(currentAngle <= 90 || currentAngle>=270)?"r":"l",
				"m",
				Color.BLACK
			);
		}
		/*
		添加图注
		*/
		int noteX = width/25;
		int noteY = height/20;
		String[] note = {
			"Criteria",
			"Criterion Ca: Reduction of Intragroup Variation",
			"Criterion Cb: Differential Metabolic Analysis",
			"Criterion Cc: Consistency in Marker Discovery",
			"Criterion Cd: Classification Accuracy"
		};
		fontSize = height/50;
		Font noteFont = new Font(null, Font.BOLD, fontSize);//用默认字体
		graphics.setFont(noteFont);
		metrics = graphics.getFontMetrics(noteFont);//字体换，标尺属性也得更新
		drawSimpleLabel(note[0], noteX, noteY, 0, 0, "l", "m", Color.BLACK);
		graphics.setFont(new Font(null, Font.BOLD, fontSize*2/3));
		for (int i=1; i < 5; i++) {
			drawSimpleBar(noteX, noteY+i*fontSize*5/4, -90, 0, fontSize, fontSize, colorReflect[4-i]);
			drawSimpleLabel(note[i], noteX + fontSize*3/2, noteY+i*fontSize*5/4, 0, 0, "l", "m", Color.BLACK);
		}
	}
	/**
	 * 绘制散射柱状图
	 * @param inputfile 输入文件名
	 * @param outputfileformat 输出文件基础名，与维度，类型一起构成完整文件名
	 * @param hasTitle 输入是否包含标题
	 * @param cutoff 筛选维度，为表现最佳的前cutoff名
	 * @param autoSize 是否自动调整尺寸，仅对jpg格式有效
	 * @param type 输出格式，暂且支持jpg/eps/pdf
	 * @throws IOException 输入输出异常
	 * @throws Exception 其他异常
	 */
	public void drawFigure(String inputfile, String outputfileformat, Boolean hasTitle, Integer cutoff, Boolean autoSize, String type) throws IOException, Exception {
		if (cutoff <=0) {
			System.err.println("cutoff should be more than zero.");
			cutoff = 100;
		}
		loadData(inputfile, hasTitle, cutoff);
		double angle = totalAngle/(criteriaList.size() - 1);
		int fontSize=56;
		int barWidth = fontSize/2;
		int r= (int) Math.round(barWidth*180/(1.5*angle*Math.PI));
		int height = 14*r;
		int width = 5*height/4;
		if (width > 23000 || !autoSize) {
			width = 23000;
			height = width*4/5;
			r = height/14;
			barWidth = (int) Math.round(Math.PI/180*angle*r*1.5);
			fontSize = barWidth*2;
		}
		type = type.toLowerCase();
		if (!type.equals("jpg")) {
			width = 2000;
			height = width*4/5;
			r = height/14;
			barWidth = (int) Math.round(Math.PI/180*angle*r*1.5);
			fontSize = barWidth*2;
		}
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
		BufferedImage image = null;
		Document document = null;
		String outfilename = String.format(outputfileformat, cutoff, type);
		switch (type) {
			case "jpg":{
				//RGBA模式不可以编码jpeg，https://community.oracle.com/message/5387869
				image = new BufferedImage(width, height, BufferedImage.TYPE_3BYTE_BGR);
				graphics = image.createGraphics();
				metrics = graphics.getFontMetrics(awtFont);
				graphics.setFont(awtFont);//设置字体
				break;
			}
			case "pdf":{
				document = new Document(new Rectangle(width, height));
				BaseFont pdfDefaultFont = BaseFont.createFont();
				FontMapper fontMapper = new FontMapper() {
					@Override
					public BaseFont awtToPdf(Font font) {
						try {
							String pdfFontName = getFontFileName(font.getFamily(), "bond");
							if (pdfFontName==null) return pdfDefaultFont;
							return BaseFont.createFont(pdfFontName, BaseFont.IDENTITY_H, BaseFont.NOT_EMBEDDED);
						} catch (Exception ex) {
							return pdfDefaultFont;
						}
					}
					@Override
					public Font pdfToAwt(BaseFont bf, int i) {
						throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
					}
				};
				PdfWriter writer = PdfWriter.getInstance(document, new FileOutputStream(outfilename));
				document.open();
				PdfContentByte pcb = writer.getDirectContent();
				graphics = pcb.createGraphics(width, height, fontMapper);
				graphics.setFont(awtFont);
				metrics = graphics.getFontMetrics(awtFont);
				
				//pcb.set
				break;
			}
			case "eps":{
				graphics = new EpsGraphics2D();
				metrics = graphics.getFontMetrics(awtFont);
				graphics.setFont(awtFont);//设置字体
				break;
			}
			default:throw new Exception("Unsupport File Type: "+ type);
		}

		graphics.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_LCD_HRGB);//文字抗锯齿
		graphics.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);//画图抗锯齿
		drawFigure(r, fontSize, barWidth, width, height, angle);
		switch(type) {
			case "jpg": {
				graphics.dispose();
				saveAsJPEG(image, outfilename, dpi);
				break;
			}
			case "pdf": {
				graphics.dispose();
				document.close();
				break;
			}
			case "eps": {
				saveAsEps((EpsGraphics2D) graphics, outfilename);
				break;
			}
			default:throw new Exception("Unsupport File Type: "+ type);
		}
	}
}
